 build_mesh 0 1 0 1 20 20
add_vertex         0.9000000        0.5000000
add_vertex         1.0000000        0.7501667
add_edge -1 -2
set_fracture -1
intersect_fractures
remove_leaves
merge_fracture_vertices .0001
polygons
mimpy_mesh mesh1.mesh
