build_mesh 0 1 0 1 20 20
add_vertex         0.0000000        0.3000000
add_vertex         1.0000000        0.3824846
add_edge -1 -2
set_fracture -1
add_vertex         0.1000000        0.3000000
add_vertex         0.5456012        0.0000000
add_edge -1 -2
set_fracture -1
add_vertex         0.1000000        0.9000000
add_vertex         0.5849643        0.6220891
add_edge -1 -2
set_fracture -1
add_vertex         0.3000000        0.1000000
add_vertex         0.0000000        0.1500462
add_edge -1 -2
set_fracture -1
add_vertex         0.5000000        0.0000000
add_vertex         0.6599994        0.4830837
add_edge -1 -2
set_fracture -1
add_vertex         0.5000000        0.4000000
add_vertex         1.0000000        0.1481836
add_edge -1 -2
set_fracture -1
add_vertex         0.8000000        0.7000000
add_vertex         0.0000000        0.9774524
add_edge -1 -2
set_fracture -1
add_vertex         1.0000000        0.1000000
add_vertex         0.1552165        0.2099880
add_edge -1 -2
set_fracture -1
intersect_fractures
remove_leaves
merge_fracture_vertices .0001
polygons
mimpy_mesh mesh1.mesh
