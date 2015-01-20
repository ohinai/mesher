
import cmd 
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib 
import pylab 
import matplotlib.patches as patches
import numpy as np
import copy 

import pickle 


import sys 

from mimpy.mesh import mesh

plt.ion()

def line_line_intersection(a1, a2, b1, b2):
    
    ts = np.linalg.lstsq(np.array([a2-a1, -b2+b1]).T, b1-a1)[0]
    
    a_intersect = a1 + (a2-a1)*ts[0]
    b_intersect = b1 + (b2-b1)*ts[1]
    
    if np.linalg.norm(a_intersect-b_intersect) > 1.e-5:
        return (None, [-1, -1])
    
    return (a_intersect, ts)

class mesher(cmd.Cmd):
    
    show_edges = True
    show_vertices = True

    show_edge_numbering = True 
    show_vertex_numbering = True
    threshold = 1.e-5

    def __init__(self): 
        cmd.Cmd.__init__(self)

        self.vertices = []
        self.edges = []
        self.vert_to_edge = {}

        self.bounding_box = [[0., 0.], [0., 0.]]

        self.polygons = None
        self.boundaries = None

        self.fracture_edges = []
        
    def emptyline(self):
        pass

    def do_toggle(self, line):
        """ Toggles what is shown in the graph.  
        """
        split_line = line.split()
        if "numbers" == split_line[0]:
            if split_line[1] == "edges":
                self.show_edge_numbering = not self.show_edge_numbering

            if split_line[1] == "vertex":
                self.show_vertex_numbering = not self.show_vertex_numbering
        
        elif split_line[0] == "edges":
            self.show_edges = not self.show_edges
            self.show_edge_numbering = False

        elif split_line[0] == "vertex":
            self.show_vertices = not self.show_vertices
            self.show_vertex_numbering = False

    def do_refresh(self, line):
        temp_axis = plt.axis()
        plt.clf()
        if self.show_vertices:
            for (index, point) in enumerate(self.vertices):
                plt.scatter(point[0], point[1])
                if self.show_vertex_numbering:
                    plt.text(point[0], point[1], str(index))

        if self.show_edges:
            for (index, edge) in enumerate(self.edges):
                point1 = self.vertices[edge[0]]
                point2 = self.vertices[edge[1]]
                plt.plot([point1[0], point2[0]], 
                         [point1[1], point2[1]], 'k-')
                if self.show_edge_numbering:
                    mid_point = (point1+point2)/2.
                    plt.text(mid_point[0], mid_point[1], str(index))
        plt.axis(temp_axis)
        plt.show()        
            
    def do_show(self, line):
        plt.clf()
        if self.show_vertices:
            for (index, point) in enumerate(self.vertices):
                plt.scatter(point[0], point[1])
                if self.show_vertex_numbering:
                    plt.text(point[0], point[1], str(index))

        if self.show_edges:
            for (index, edge) in enumerate(self.edges):
                point1 = self.vertices[edge[0]]
                point2 = self.vertices[edge[1]]
                plt.plot([point1[0], point2[0]], 
                         [point1[1], point2[1]], 'k-')
                if self.show_edge_numbering:
                    mid_point = (point1+point2)/2.
                    plt.text(mid_point[0], mid_point[1], str(index))
        plt.show()        
        

    def find_orphaned_edges(self ):
        """ Finds all edges are either not connected 
        at all, or connected only on one side. 
        """
        single_sided = set()
        no_sided = set()
        for vertex_index in range(len(self.vertices)):
            if self.vert_to_edge.has_key(vertex_index):
                if len(self.vert_to_edge[vertex_index]) ==  1:
                    edge_index = self.vert_to_edge[vertex_index][0]
                    if edge_index in single_sided:
                        single_sided.remove(edge_index)
                        no_sided.add(edge_index)
                    else:
                        single_sided.add(edge_index)
        return [single_sided, no_sided]

    def do_find_orphans(self, line):
        """ Finds all edges that either unconnected to 
        anything or conntected only on a single side. 
        """
        [single_sided, no_sided] = self.find_orphaned_edges()
        print "single sided edges", single_sided
        print "no sided edges", no_sided
        
    def find_next_edge_cc(self, edge_index, vertex):
        """ Finds next edge looping 
        couter-clockwise. 
        """
        local_index = self.vert_to_edge[vertex].index(edge_index)
        next_edge = self.vert_to_edge[vertex][(local_index+1)%len(self.vert_to_edge[vertex])]
        return next_edge

    def find_next_edge_c(self, edge_index, vertex):
        """ Finds next edge looping 
        couter-clockwise. 
        """
        local_index = self.vert_to_edge[vertex].index(edge_index)
        next_edge = self.vert_to_edge[vertex][local_index-1]
        return next_edge

    def do_polygons(self, line):
        """ Compute all polgons in the graph.
        """
        done_edges = set()
        polygons = []
        boundaries = []
        paths = []
        index_paths = []
        def direction_int(v1, v2):
            point1 = self.vertices[v1]
            point2 = self.vertices[v2]
            return (point2[0]-point1[0])*(point2[1]+point1[1])

        for edge_index in range(len(self.edges)):
            if (edge_index, 1) not in done_edges:
                direction_sum = 0.
                current_polygon = [edge_index]
                [current_v1, current_v2] = self.edges[edge_index]
                original_v1 = current_v1
                current_path = [self.vertices[current_v1], self.vertices[current_v2]]
                current_index_path = [(edge_index, 1)]
                direction_sum += direction_int(current_v1, current_v2)
                ## Loop counter-clockwise
                next_edge = self.find_next_edge_cc(edge_index, current_v2)

                current_polygon.append(next_edge)
                [next_v1, next_v2] = self.edges[next_edge]

                if next_v1 == current_v2:
                    next_vertex = next_v2
                    direction_sum += direction_int(next_v1, next_v2)
                    done_edges.add((next_edge, 1))
                    current_index_path.append((next_edge, 1))

                elif next_v2 == current_v2:
                    next_vertex = next_v1
                    direction_sum += direction_int(next_v2, next_v1)
                    done_edges.add((next_edge, -1))
                    current_index_path.append((next_edge, -1))

                else:
                    raise Exception("Problem traversing graph")
                while next_vertex != original_v1:
                    current_path.append(self.vertices[next_vertex])

                    next_edge = self.find_next_edge_cc(next_edge, next_vertex)

                    current_polygon.append(next_edge)
                    [next_v1, next_v2] = self.edges[next_edge]

                    if next_v1 == next_vertex:
                        next_vertex = next_v2
                        direction_sum += direction_int(next_v1, next_v2)
                        done_edges.add((next_edge, 1))
                        current_index_path.append((next_edge, 1))

                    elif next_v2 == next_vertex:
                        next_vertex = next_v1
                        direction_sum += direction_int(next_v2, next_v1)
                        done_edges.add((next_edge, -1))
                        current_index_path.append((next_edge, -1))

                    else:
                        raise Exception("Problem traversing graph")

                if direction_sum < 0:
                    boundaries += current_polygon
                else:
                    polygons.append(current_polygon)
                    paths.append(current_path)
                    index_paths.append(current_index_path)

            if (edge_index, -1) not in done_edges:
                ## Loop clockwise                
                direction_sum = 0.

                current_polygon = [edge_index]
                [current_v1, current_v2] = self.edges[edge_index]
                original_v1 = current_v1
                current_path = [self.vertices[current_v1], self.vertices[current_v2]]
                current_index_path = [(edge_index, -1)]
                direction_sum += direction_int(current_v1, current_v2)
                next_edge = self.find_next_edge_c(edge_index, current_v2)

                current_polygon.append(next_edge)
                [next_v1, next_v2] = self.edges[next_edge]
                if next_v1 == current_v2:
                    next_vertex = next_v2
                    direction_sum += direction_int(next_v1, next_v2)
                    done_edges.add((next_edge, -1))
                    current_index_path.append((next_edge, -1))

                elif next_v2 == current_v2:
                    next_vertex = next_v1
                    direction_sum += direction_int(next_v2, next_v1)
                    done_edges.add((next_edge, 1))
                    current_index_path.append((next_edge, 1))

                else:
                    raise Exception("Problem traversing graph")
                while next_vertex != original_v1:
                    current_path.append(self.vertices[next_vertex])
                    
                    next_edge = self.find_next_edge_c(next_edge, next_vertex)

                    current_polygon.append(next_edge)
                    [next_v1, next_v2] = self.edges[next_edge]
                    if next_v1 == next_vertex:
                        next_vertex = next_v2
                        direction_sum += direction_int(next_v1, next_v2)
                        done_edges.add((next_edge, -1))
                        current_index_path.append((next_edge, -1))
                    elif next_v2 == next_vertex:
                        next_vertex = next_v1                        
                        direction_sum += direction_int(next_v2, next_v1)
                        done_edges.add((next_edge, 1))
                        current_index_path.append((next_edge, 1))

                    else:
                        raise Exception("Problem traversing graph")

                if direction_sum > 0:
                    boundaries += current_polygon
                else:
                    polygons.append(current_polygon)
                    paths.append(current_path)
                    current_index_path.reverse()
                    index_paths.append(current_index_path)
                    
        self.polygons = index_paths
        self.boundaries = boundaries

        fig=pylab.figure()
        ax=fig.add_subplot(111)
        patches = map(lambda x:Polygon(x, True), map(lambda x:np.array(x), paths))
        p = PatchCollection(patches,cmap=matplotlib.cm.jet,  alpha=1.0)

#        p.set_array(np.random.rand(len(paths)))
        p.set_array(np.ones(len(paths)))
        ax.add_collection(p)
        p.set_clim([0., 1.])

        ax.autoscale_view(True, True, True)

        pylab.show()
    
    def do_set_fracture(self, line):
        """ Sets edge as fracture edge. 
        """
        edge_index = int(line.split()[0])
        self.fracture_edges.append(edge_index)
        
    def do_mimpy_mesh(self, line):
        """ Produces a mimpy mesh from the 
        graph and stores in the filename specified.  
        """
        res_mesh = mesh.Mesh()
        
        file_name = line
        res_mesh.dim = 3

        for point in self.vertices:
            res_mesh.add_point(np.array([point[0], point[1], 0.]))
        
        for (index, point) in enumerate(self.vertices):
            new_index = res_mesh.add_point(np.array([point[0],point[1], -1.]))

        done_edges = set()
        edge_to_face_map = {}
        
        jump = len(self.vertices)

        res_mesh.add_boundary_marker(0, "bottom surface")
        res_mesh.add_boundary_marker(1, "top surface")
        
        res_mesh.add_boundary_marker(2, "perimeter")
        
        for polygon in self.polygons:
            current_cell = []
            current_cell_orientations = []
            top_face = []
            bot_face = []
            for (edge_index, direction) in polygon:
                if direction > 0:
                    [p1, p2]= self.edges[edge_index]
                else:
                    [p2, p1]= self.edges[edge_index]
                    
                top_face.append(p1)
                bot_face.append(p1+jump)
               
                if edge_index in done_edges:
                    current_cell.append(edge_to_face_map[edge_index])
                    current_cell_orientations.append(-1)
                    
                else:
                    new_face = [p1, p2, p2+jump, p1+jump]
                    new_face_index = res_mesh.add_face(new_face)

                    normal = res_mesh.find_face_normal(new_face_index)
                    res_mesh.set_face_normal(new_face_index, normal)
                    (area, centroid) = res_mesh.find_face_centroid(new_face_index)
                    res_mesh.set_face_area(new_face_index, area)
                    res_mesh.set_face_real_centroid(new_face_index, centroid)

                    current_cell.append(new_face_index)
                    current_cell_orientations.append(1)
                    edge_to_face_map[edge_index] = new_face_index
                    done_edges.add(edge_index)

            top_face_index = res_mesh.add_face(top_face)

            normal = res_mesh.find_face_normal(top_face_index)
            res_mesh.set_face_normal(top_face_index, normal)
            (area, centroid) = res_mesh.find_face_centroid(top_face_index)
            res_mesh.set_face_area(top_face_index, area)
            res_mesh.set_face_real_centroid(top_face_index, centroid)

            res_mesh.set_face_quadrature_points(top_face_index, 
                                                [centroid])
            res_mesh.set_face_quadrature_weights(top_face_index, 
                                                 [area])

            bot_face_index = res_mesh.add_face(bot_face)
            
            normal = res_mesh.find_face_normal(bot_face_index)
            res_mesh.set_face_normal(bot_face_index, normal)
            (area, centroid) = res_mesh.find_face_centroid(bot_face_index)
            res_mesh.set_face_area(bot_face_index, area)
            res_mesh.set_face_real_centroid(bot_face_index, centroid)

            res_mesh.set_face_quadrature_points(bot_face_index, 
                                                [centroid])
            res_mesh.set_face_quadrature_weights(bot_face_index, 
                                                 [area])

            res_mesh.add_boundary_face(0,  bot_face_index, -1)
            res_mesh.add_boundary_face(1,  top_face_index, 1)

            current_cell.append(top_face_index)
            current_cell_orientations.append(1)
            current_cell.append(bot_face_index)
            current_cell_orientations.append(-1)
            
            new_cell_index = res_mesh.add_cell(current_cell, current_cell_orientations)
            
            res_mesh.set_cell_k(new_cell_index, 1.e-4*np.eye(3))
            
            (volume, centroid) = res_mesh.find_volume_centroid(new_cell_index)
            res_mesh.set_cell_real_centroid(new_cell_index, centroid)
            res_mesh.set_cell_volume(new_cell_index, volume)
            

        for edge_index in self.boundaries:
            face_index = edge_to_face_map[edge_index]
            res_mesh.add_boundary_face(2, face_index, 1)
            res_mesh.set_face_quadrature_points(face_index, 
                                                [res_mesh.get_face_real_centroid(face_index)])
            res_mesh.set_face_quadrature_weights(face_index, 
                                                 [res_mesh.get_face_area(face_index)])
            

        res_mesh.build_frac_from_faces([edge_to_face_map[edge_index] for edge_index in self.fracture_edges])
            
        res_mesh.output_vtk_mesh(file_name, [res_mesh.get_cell_domain_all(), 
                                             [k[0, 0] for k in res_mesh.get_all_k()]], ["DOMAIN", "K"])
        pickle_file = open(file_name, 'w')
        pickle.dump(res_mesh, pickle_file)
        pickle_file.close()

    def do_EOF(self, line):
        return True

    def add_vertex(self, point):
        """ Adds new vertex, returns 
        vertex index. 
        """
        self.vertices.append(point)
        self.bounding_box[0][0] = min(self.bounding_box[0][0], point[0])
        self.bounding_box[1][0] = max(self.bounding_box[0][0], point[0])

        self.bounding_box[0][1] = min(self.bounding_box[0][1], point[1])
        self.bounding_box[1][1] = max(self.bounding_box[0][1], point[1])

        return len(self.vertices)-1

    def do_add_vertex(self, line):
        new_point = np.array(map(float, line.split()))
        new_index = self.add_vertex(new_point)
        print "new vertex", "[", new_index, "]", new_point

    def find_cc_position(self, index_v1, index_v2, edges):
        """ Given edge p1->p2, find position so that 
        all edges around p2 are counter-clockwise order. 
        """
        points = []
        if len(edges)<1:
            return 0
        for edge_index in edges:
            [current_v1, current_v2] = self.edges[edge_index]
            if current_v1 == index_v2:
                points.append(self.vertices[current_v2])
            elif current_v2 == index_v2:
                points.append(self.vertices[current_v1])
            else:
                print "vert_to_edge error"

        p1 = self.vertices[index_v1]
        p2 = self.vertices[index_v2]

        signs = []
        v1 = p2-p1
        all_pos = True
        all_neg = True
        for point in points:
            v2 = p2-point
            if np.cross(np.array([v1[0], v1[1], 0.]),
                        np.array([v2[0], v2[1], 0.]))[2] > 0.:
                signs.append(1)
                all_neg = False
            else:
                signs.append(-1)
                all_pos = False

        if all_pos:
            return 0
        if all_neg:
            return len(points) 
        if signs[-1]<0 and signs[0]>0:
            return 0 
        found_negative = False
        for (index, sign) in enumerate(signs):
            if sign < 0:
                found_negative = True
            if sign > 0 and found_negative:
                return index

    def do_show_v2e(self, line):
        """ Shows the edges connected with 
        a vertex.  
        """
        plt.clf()
        v_index = int(line.split()[0])
        print "showing vertices connected with", v_index
        vertex = self.vertices[v_index]
        plt.scatter(vertex[0], vertex[1])
        for (loc, edge_index) in enumerate(self.vert_to_edge[v_index]):
            [v1, v2] = self.edges[edge_index]
            p1 = self.vertices[v1]
            p2 = self.vertices[v2]
            plt.plot([p1[0], p2[0]], [p1[1], p2[1]], 'k-')
            center = (p2+p1)/2.
            plt.text(center[0], center[1], str(loc))
        plt.show()
        
    def add_edge(self, v1, v2):
        self.edges.append([v1, v2])

        self.update_v_to_e(len(self.edges)-1, v1, v2)
        return len(self.edges)-1


    def do_add_edge(self, line):
        [v1, v2] = map(int, line.split())
        edge_index = self.add_edge(v1, v2)
        print "new edge","[", edge_index , "]", v1, "<->", v2

    def update_v_to_e(self, edge_index, v1, v2 ):
        """ Updates the vertex to edge map 
        for new edge. 
        """
        if self.vert_to_edge.has_key(v1):
            pos = self.find_cc_position(v2, v1, self.vert_to_edge[v1])
            self.vert_to_edge[v1].insert(pos, edge_index)
        else:
            self.vert_to_edge[v1] = [len(self.edges)-1]

        if self.vert_to_edge.has_key(v2):
            pos = self.find_cc_position(v1, v2, self.vert_to_edge[v2])
            self.vert_to_edge[v2].insert(pos, edge_index)
        else:
            self.vert_to_edge[v2] = [edge_index]
            

    def set_edge(self, edge_index, v1, v2):
        """ Modifies vertices of existing 
        edge. 
        """
        [orig_v1, orig_v2] = self.edges[edge_index]

        self.vert_to_edge[orig_v1].remove(edge_index)
        self.vert_to_edge[orig_v2].remove(edge_index)

        self.edges[edge_index][0] = v1
        self.edges[edge_index][1] = v2

        self.update_v_to_e(edge_index, v1, v2)

    def intersect_edge(self, edge_index):
        """ For a given edge number, intersect it 
        with all other edges in the graph. 
        """
        edge = self.edges[edge_index]
        point1 = self.vertices[edge[0]]
        point2 = self.vertices[edge[1]]
        
        # edges in graph to be modified 
        edge_mod = []
        
        # vertices on the primary line
        segments = []
        for (current_index, current_edge) in enumerate(self.edges):
            if edge_index == current_index:
                pass
            else:
                current_point1 = self.vertices[current_edge[0]]
                current_point2 = self.vertices[current_edge[1]]
                
                (intersection, param) = line_line_intersection(point1,
                                                               point2,
                                                               current_point1,
                                                               current_point2)
                if (abs(param[0])<self.threshold or abs(param[0]-1)<self.threshold) and \
                        (abs(param[1])<self.threshold or abs(param[1]-1)<self.threshold):
                    pass
                elif 0-self.threshold<=param[0]<=1.+self.threshold \
                        and 0.-self.threshold<=param[1]<=1.+self.threshold:
                    if abs(param[0])<self.threshold:
                        edge_mod.append([current_index, edge[0]])
                    elif abs(param[0]-1)<self.threshold:
                        edge_mod.append([current_index, edge[1]])

                    else:
                        if abs(param[1]) < self.threshold:
                            segments.append([param[0], current_edge[0]])
                        elif  abs(param[1]-1.) < self.threshold:
                            segments.append([param[0], current_edge[1]])
                        else:
                            new_point = self.add_vertex(intersection)
                            edge_mod.append([current_index, new_point])
                            segments.append([param[0], new_point])

        for current_mod in edge_mod:
            current_index = current_mod[0]
            point_index = current_mod[1]
            self.add_edge(self.edges[current_index][1], point_index)
            self.set_edge(current_index, self.edges[current_index][0], point_index )
            
        segments.sort()

        segments.append([1. , self.edges[edge_index][1]])
        self.set_edge(edge_index, self.edges[edge_index][0], segments[0][1])
        done_edges = set()
        done_edges.add(edge_index)

        for seg_index in range(len(segments)-1):
            v1 = segments[seg_index][1]
            v2 = segments[seg_index+1][1]
            new_edge_index = self.add_edge(v1, v2)
            done_edges.add(new_edge_index)
        return done_edges

    def update_edge_numbering(self, edge_index):
        """ After removing an edge, updates 
        the edge numbering everywhere that is needed. 
        """
        for vertex_index in self.vert_to_edge:
            for (index,  current_edge_index) in enumerate(self.vert_to_edge[vertex_index]):
                if current_edge_index > edge_index:
                    self.vert_to_edge[vertex_index][index] -= 1

    def do_remove_edge(self, line):
        """ Removes edge from graph. 
        """
        index = int(line.split()[0])
        [v1, v2] = self.edges[index]
        self.vert_to_edge[v1].remove(index)
        self.vert_to_edge[v2].remove(index)
        self.edges.pop(index)
        self.update_edge_numbering(index)
        
    def do_intersect_all(self, line):
        """ Intersect all edges. 
        """

        current_index = 0
        done_edges = set() 
        while current_index < len(self.edges):
            
            if current_index not in done_edges:
                done_edges.union(self.intersect_edge(current_index))

            current_index += 1
        
        
    def do_intersect_edge(self, line):
        """ For a given edge number, intersect it 
        with all other edges in the graph. 
        """
        line_split = line.split()
        edge_index = int(line_split[0])
        self.intersect_edge(edge_index)

    def do_load_commands(self, line):
        file = open(line)
        for line in file:
            if len(line.split())>0:
                self.onecmd(line)

    def do_build_mesh(self, line):
        """ Builds rectangular mesh. Input:
        [x dimension]
        [y dimension]
        [number of cells in x]
        [number of cells in y]
        """
        split_line = line.split()
        x_dim = float(split_line[0])
        y_dim = float(split_line[1])
        
        nx = int(split_line[2])
        ny = int(split_line[3])
        ij_to_point = {}

        dx = x_dim/nx
        dy = y_dim/ny

        for i in range(nx+1):
            for j in range(ny+1):
                new_point = np.array([i*dx, j*dy])
                new_point_index = self.add_vertex(new_point)
                ij_to_point[(i, j)] = new_point_index

        for i in range(nx+1):
            for j in range(ny+1):
                
                if i < nx:
                    self.add_edge(ij_to_point[(i, j)], 
                                  ij_to_point[(i+1, j)])
                
                if j < ny:
                    self.add_edge(ij_to_point[(i, j)], 
                                  ij_to_point[(i, j+1)])
                    
        
mesher().cmdloop()
