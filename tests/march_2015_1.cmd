build_mesh 0 1 0 1 20 20
add_vertex         0.2000000        0.2000000
add_vertex         0.6070302        0.8076009
add_edge -1 -2
set_fracture -1
add_vertex         0.8000000        0.4000000
add_vertex         0.4133587        0.6610174
add_edge -1 -2
set_fracture -1
add_vertex         0.8000000        0.8000000
add_vertex         0.8850909        1.0000000
add_edge -1 -2
set_fracture -1
add_vertex         0.8000000        0.9000000
add_vertex         0.5445319        1.0000000
add_edge -1 -2
set_fracture -1
add_vertex         0.9000000        0.7000000
add_vertex         0.9394587        0.9653175
add_edge -1 -2
set_fracture -1
add_vertex         0.9000000        1.0000000
add_vertex         0.0000000        0.9993066
add_edge -1 -2
set_fracture -1
intersect_fractures
remove_leaves
merge_fracture_vertices .0001
polygons
mimpy_mesh mesh1.mesh

