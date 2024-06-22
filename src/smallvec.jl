# Inline a small vector into a struct. Meaning use
# a StaticVector for the first N elements and a
# Vector allocated as needed for the rest. This
# is to reduce allocations because many facets have
# very few 'above' points.

using MacroTools, StaticArrays

# todo: make SmallVec work with type params from struct header?
macro embed_SmallVec(expr)
    def = MacroTools.splitstructdef(expr)
    !def[:mutable] && error("Expected mutable struct")

    vectors = []
    seen_sv_field = false
    fields′ = map(def[:fields]) do (fieldname, fieldtype)
        match = @capture(fieldtype, @SmallVec(params_))
        if !match
            seen_sv_field && error("SmallVec fields must be last")
            return [(fieldname, fieldtype)]
        end
        seen_sv_field = true

        match = @capture(params, {eltype_, rest__})
        if !match || length(rest) > 1 || (length(rest) == 1 && !(rest[1] isa Int))
            error("Expected @SmallVec{type} or @SmallVec{type, N}")
        end

        push!(vectors, (fieldname, eltype))
        N = isempty(rest) ? 4 : rest[1]
        return [
            (Symbol(fieldname, :_head),   :(SVector{$N, $eltype}))
            (Symbol(fieldname, :_tail),   :(Union{Nothing, Vector{$eltype}}))
            (Symbol(fieldname, :_length), :Int)
        ]
    end
    def[:fields] = reduce(vcat, fields′)

    view_structs = []
    view_accessors = []
    initializers = []
    for (name, eltype) in vectors
        view_name = Symbol(def[:name], :_, name, :_view)

        fhead = Symbol(name, :_head)
        ftail = Symbol(name, :_tail)
        flength = Symbol(name, :_length)

        accessor = quote
            (sym === $(QuoteNode(name))) && return $view_name(x)
        end
        push!(view_accessors, accessor)

        ex = quote
            struct $view_name{T} <: AbstractVector{$eltype}
                obj::T
            end

            function Base.iterate(view::$view_name, state=nothing)
                idx = isnothing(state) ? 1 : state + 1

                if idx <= view.obj.$flength
                    if idx <= length(view.obj.$fhead)
                        return (@inbounds view.obj.$fhead[idx]), idx
                    else
                        tail = view.obj.$ftail
                        if isnothing(tail)
                            return nothing
                        else
                            (@inbounds tail[idx - length(view.obj.$fhead)]), idx
                        end
                    end
                else
                    return nothing
                end
            end

            Base.axes(view::$view_name) = (Base.OneTo(view.obj.$flength), )
            Base.size(view::$view_name) = (view.obj.$flength, )
            Base.eltype(view::$view_name) = $eltype

            function Base.push!(view::$view_name, el)
                if view.obj.$flength < length(view.obj.$fhead)
                    view.obj.$fhead = setindex(view.obj.$fhead, el, view.obj.$flength + 1)
                else
                    if isnothing(view.obj.$ftail)
                        view.obj.$ftail = [el]
                    else
                        push!(view.obj.$ftail, el)
                    end
                end
                view.obj.$flength += 1
                return view
            end

            function Base.setindex!(view::$view_name, el, i)
                #@boundscheck checkbounds(Bool, view, i)

                if i <= length(view.obj.$fhead)
                    view.obj.$fhead = @inbounds setindex(view.obj.$fhead, el, i)
                else
                    @inbounds view.obj.$ftail[i - length(view.obj.$fhead)] = el
                end
                return view
            end

            function Base.getindex(view::$view_name, i)
                if i <= length(view.obj.$fhead)
                    return @inbounds view.obj.$fhead[i]
                else
                    return @inbounds view.obj.$ftail[i - length(view.obj.$fhead)]
                end
            end

            function Base.empty!(view::$view_name)
                view.obj.$flength = 0
                if !isnothing(view.obj.$ftail)
                    empty!(view.obj.$ftail)
                end
                return view
            end
        end
        push!(view_structs, ex)

        init = quote
            __self.$fhead = zero(fieldtype(typeof(__self), $(QuoteNode(fhead))))
            __self.$ftail = nothing
            __self.$flength = 0
        end
        push!(initializers, init)
    end

    def[:constructors] = map(def[:constructors]) do constructor
        MacroTools.postwalk(constructor) do x
            @capture(x, new(args__)) || @capture(x, new{ts__}(args__)) || return x
            
            return quote
                __self = $x
                $(Expr(:block, initializers...))
                __self
            end
        end
    end

    getproperty_fn = quote
        function Base.getproperty(x::$(def[:name]), sym::Symbol)
            $(Expr(:block, view_accessors...))

            return Base.getfield(x, sym)
        end
    end

    res = quote
        $(MacroTools.combinestructdef(def))
        $getproperty_fn

        $(Expr(:block, view_structs...))
    end

    return esc(res)
end
