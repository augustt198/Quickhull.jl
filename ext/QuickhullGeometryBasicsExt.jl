module QuickhullGeometryBasicsExt

using Quickhull, GeometryBasics
using MappedArrays

function GeometryBasics.Mesh(hull::Quickhull.Hull{D, T, K, V}) where {D, T, K, V}
    P = Point{D, T}
    F = NgonFace{D, Int}

    points = mappedarray(hull.pts) do point
        Point{D, T}(point)
    end
    
    faces = mappedarray(hull.facets.arr) do facet
        NgonFace{D, Int}(facet.plane.point_indices)
    end

    fv = FaceView{GeometryBasics.Ngon{D, T, D, P}, P, F, typeof(points), typeof(faces)}(points, faces)
    GeometryBasics.Mesh(fv)
end

end # module
