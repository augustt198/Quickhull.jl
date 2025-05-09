import PrecompileTools

PrecompileTools.@setup_workload begin
    pts_vec2d = [Point2(sincos(θ)...) for θ = 0:0.1:2π]
    pts_mat2d = reduce(hcat, pts_vec2d)

    pts_vec3d = [Point3(sincos(θ)..., θ) for θ = 0:0.1:2π]
    pts_mat3d = reduce(hcat, pts_vec3d)

    PrecompileTools.@compile_workload begin
        quickhull(pts_vec2d)
        quickhull(pts_mat2d)
        quickhull(pts_vec3d)
        quickhull(pts_mat3d)

        delaunay(pts_vec2d)
        delaunay(pts_mat2d)
        delaunay(pts_vec3d)
        delaunay(pts_mat3d)
    end

end
