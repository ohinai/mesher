#! /usr/bin/env python
import cmd 
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.font_manager import FontProperties
import matplotlib 
import pylab 
import matplotlib.patches as patches
import numpy as np
import copy 

import pickle 
import matplotlib.image as mpimg 

import sys 
import io

import random

from mimpy.mesh import mesh

plt.ion()

def line_line_intersection(a1, a2, b1, b2):
    
    ts = np.linalg.lstsq(np.array([a2-a1, -b2+b1]).T, b1-a1)[0]
    
    a_intersect = a1 + (a2-a1)*ts[0]
    b_intersect = b1 + (b2-b1)*ts[1]
    
    if np.linalg.norm(a_intersect-b_intersect) > 1.e-5:
        return (None, [-1, -1])
    
    return (a_intersect, ts)

def is_point_on_ray(point, a1, a2):
    v1 = a2-a1
    v2 = point-a1
    if np.linalg.norm(v2)<1.e-8:
        return True
    if np.linalg.norm(point-a2)<1.e-8:
        return True
    
    v2_dot_v1 = np.linalg.norm(v2.dot(v1/np.linalg.norm(v1)))
    if abs(v2_dot_v1-np.linalg.norm(v2))<1.e-8:
        return True
    return False
    
def is_point_on_line(point, a1, a2):
    v1 = a2-a1
    v2 = point-a1
    if np.linalg.norm(v1)<np.linalg.norm(v2):
        return False
    if np.linalg.norm(v2)<1.e-8:
        return True
    if np.linalg.norm(point-a2)<1.e-8:
        return True
    v2_dot_v1 = np.linalg.norm(v2.dot(v1/np.linalg.norm(v1)))

    if abs(v2_dot_v1-np.linalg.norm(v2))<1.e-8:
        return True
    return False

class Mesher(cmd.Cmd):
    
    show_edges = True
    show_vertices = True

    show_edge_numbering = True 
    show_vertex_numbering = True

    def __init__(self): 
        cmd.Cmd.__init__(self)

        self.threshold = 1.e-5
        self.point_threshold = 1.e-3
        
        self.vertices = []
        self.edges = []
        self.vert_to_edge = {}

        self.bounding_box = [[0., 0.], [0., 0.]]

        self.holes = []

        self.polygons = None
        self.boundaries = None

        self.fracture_edges = []

        self.session = []
        
        self.highlight_edges = None

        self.reservoir_depth = 1.

        ## Maintains the indices of the last 
        ## vertices to be added to the graph. 
        ## Used for negative indexing. 
        self.vertex_cache = []
        self.vertex_cache_size = 10

    def emptyline(self):
        pass

    def precmd(self, line):
        self.session.append(line)
        if len(line)>0 and line[0] == "#":
            return ""
        return line

    def ee_intersection(self, edge_1, edge_2):
        """ Finds intersection of two edges. 
        Returns interesction point and tuple 
        representing the parametric value 
        of point on teach line. 
        """
        v1 = self.vertices[self.edges[edge_1][0]]
        v2 = self.vertices[self.edges[edge_1][1]]

        v3 = self.vertices[self.edges[edge_2][0]]
        v4 = self.vertices[self.edges[edge_2][1]]

        return line_line_intersection(v1, v2, v3, v4)

    def do_save_mesh(self, line):
        """ Saves current mesh in file. 
        """
        file_name = line.split()[0]
        
        output_file = open(file_name, 'w')
        print >>output_file, "VERTICES", len(self.vertices)
        for vertex in self.vertices:
            print >>output_file, vertex[0], vertex[1]

        print >>output_file, "EDGES", len(self.edges)
        for edge in self.edges:
            print >>output_file, edge[0], edge[1]

        print >>output_file, "VTOE", len(self.vert_to_edge)
        for key in self.vert_to_edge:
            print >>output_file, key, " ".join(map(str, self.vert_to_edge[key]))
        
        print >>output_file, "FRACTURE", len(self.fracture_edges)
        for edge in self.fracture_edges:
            print >>output_file, edge
        
    def do_load_mesh(self, line):
        """ Loads mesh from file. 
        """
        file_name = line.split()[0]
        
        output_file = open(file_name)
        
        current_line = output_file.readline()
        line_split = current_line.split()
        n_vertices = int(line_split[1])
        for v_index in range(n_vertices):
            current_line = output_file.readline()
            line_split = current_line.split()
            current_point = np.array(map(float, line_split))
            self.vertices.append(current_point)            

        current_line = output_file.readline()
        line_split = current_line.split()
        n_edges = int(line_split[1])
        
        for e_index in range(n_edges):
            current_line = output_file.readline()
            line_split = current_line.split()
            current_edge = np.array(map(int, line_split))
            self.edges.append(current_edge)

        current_line = output_file.readline()
        line_split = current_line.split()
        n_vtoe = int(line_split[1])
        
        for index in range(n_vtoe):
            current_line = output_file.readline()
            line_split = current_line.split()

            self.vert_to_edge[int(line_split[0])] = map(int, line_split[1:])
            
 
        current_line = output_file.readline()
        line_split = current_line.split()
        n_fractures = int(line_split[1])
        
        for index in range(n_fractures):
            current_line = output_file.readline()
            line_split = current_line.split()
            
            self.fracture_edges.append(int(line_split[0]))
            

    def do_save_session(self, line):
        """ Save session in file. 
        """
        output_file = open(line,  'w')
        for line in self.session:
            print>>output_file, line
        output_file.close()
        
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
            for edge_index in self.fracture_edges:
                edge = self.edges[edge_index]
                point1 = self.vertices[edge[0]]
                point2 = self.vertices[edge[1]]
                plt.plot([point1[0], point2[0]], 
                         [point1[1], point2[1]], 'k-', color = 'r', linewidth = 7.)

            for (index, edge) in enumerate(self.edges):
                point1 = self.vertices[edge[0]]
                point2 = self.vertices[edge[1]]
                plt.plot([point1[0], point2[0]], 
                         [point1[1], point2[1]], 'k-')
                if self.show_edge_numbering:
                    mid_point = (point1+point2)/2.
                    plt.text(mid_point[0], mid_point[1], str(index))

            if self.highlight_edges != None:
                edge = self.edges[self.highlight_edges]
                point1 = self.vertices[edge[0]]
                point2 = self.vertices[edge[1]]
                plt.plot([point1[0], point2[0]], 
                         [point1[1], point2[1]], 'k-', color = 'y', linewidth = 8.)

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
            for edge_index in self.fracture_edges:
                edge = self.edges[edge_index]
                point1 = self.vertices[edge[0]]
                point2 = self.vertices[edge[1]]
                plt.plot([point1[0], point2[0]], 
                         [point1[1], point2[1]], 'k-', color = 'r', linewidth = 7.)
                

            for (index, edge) in enumerate(self.edges):
                point1 = self.vertices[edge[0]]
                point2 = self.vertices[edge[1]]
                plt.plot([point1[0], point2[0]], 
                         [point1[1], point2[1]], 'k-')
                if self.show_edge_numbering:
                    mid_point = (point1+point2)/2.
                    plt.text(mid_point[0], mid_point[1], str(index))

            if self.highlight_edges != None:
                edge = self.edges[self.highlight_edges]
                point1 = self.vertices[edge[0]]
                point2 = self.vertices[edge[1]]
                plt.plot([point1[0], point2[0]], 
                         [point1[1], point2[1]], 'k-', color = 'y', linewidth = 8.)

        plt.axis('off')
        plt.show()

    def do_save_pdf(self, line):
        """ Saves the figure as a pdf. 
        """
        plt.clf()
        font1 = FontProperties()
        font1.set_size(1)
        if self.show_vertices:
            for (index, point) in enumerate(self.vertices):
                plt.scatter(point[0], point[1])
                if self.show_vertex_numbering:
                    plt.text(point[0], point[1], str(index), fontproperties=font1)

        if self.show_edges:
            for edge_index in self.fracture_edges:
                edge = self.edges[edge_index]
                point1 = self.vertices[edge[0]]
                point2 = self.vertices[edge[1]]
                plt.plot([point1[0], point2[0]], 
                         [point1[1], point2[1]], 'k-', color = 'r', linewidth = 7.)
                
            for (index, edge) in enumerate(self.edges):
                point1 = self.vertices[edge[0]]
                point2 = self.vertices[edge[1]]
                plt.plot([point1[0], point2[0]], 
                         [point1[1], point2[1]], 'k-')
                if self.show_edge_numbering:
                    mid_point = (point1+point2)/2.
                    plt.text(mid_point[0], mid_point[1], str(index))

            if self.highlight_edges != None:
                edge = self.edges[self.highlight_edges]
                point1 = self.vertices[edge[0]]
                point2 = self.vertices[edge[1]]
                plt.plot([point1[0], point2[0]], 
                         [point1[1], point2[1]], 'k-', color = 'y', linewidth = 8.)

        plt.axis('off')
        name = line.split()[0]
        pylab.savefig(name+".pdf")
        
    def do_remove_leaves(self, line):
        """ Remove all leaf edges. 
        """
        (single_sided, no_sided) = self.find_leaf_edges()
        to_be_removed = list(single_sided) + list(no_sided)
        to_be_removed.sort()
        to_be_removed.reverse()
        for edge_index in to_be_removed:
            self.remove_edge(edge_index)

    def do_mo_evtoe(self, line):
        """ Shifts vertex edge to by on existing edge. 
        """
        line_split = line.split()
        edge_index = int(line_split[0])
        vertex_index = int(line_split[1])

        target_edge = int(line_split[2])
        
        original = self.vertices[vertex_index]
        (point, ts) = self.ee_intersection(edge_index, target_edge)
        
        self.vertices[vertex_index] = point
        print "vertex", vertex_index, "shifted:", original, "->", point

    def do_ex_evtoe(self, line):
        """ Extends edge to another edge by 
        finding the inersection on the second edge 
        and adding a new edge between. 
        """
        line_split = line.split()
        edge_index = int(line_split[0])
        vertex_index = int(line_split[1])

        target_edge = int(line_split[2])
        
        (point, ts) = self.ee_intersection(edge_index, target_edge)
        
        new_vertex=self.add_vertex(point)
        new_edge = self.add_edge(vertex_index, new_vertex)
        print "new vertex", "[", new_vertex, "]", point 
        print "new edge", "[", new_edge, "]",  vertex_index, "<->", new_vertex


    def find_leaf_edges(self ):
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

    def do_find_leaves(self, line):
        """ Finds all edges that either unconnected to 
        anything or conntected only on a single side. 
        """
        [single_sided, no_sided] = self.find_leaf_edges()
        print "single sided edges", single_sided
        print "no sided edges", no_sided
        
    def find_next_edge_cc(self, edge_index, vertex):
        """ Finds next edge looping 
        couter-clockwise. 
        """
        local_index = self.vert_to_edge[vertex].index(edge_index)
        next_edge_index = self.vert_to_edge[vertex][(local_index+1)%len(self.vert_to_edge[vertex])]
        return next_edge_index

    def find_next_edge_vertex_cc(self, edge_index, vertex):
        """ Finds next edge and vertex looping 
        couter-clockwise. 
        """
        #local_index = self.vert_to_edge[vertex].index(edge_index)
        #next_edge_index = self.vert_to_edge[vertex][(local_index+1)%len(self.vert_to_edge[vertex])]
        #next_edge = self.edges[next_edge_index]
        local_index = self.vert_to_edge[vertex].index(edge_index)
        next_edge_index = self.vert_to_edge[vertex][local_index-1]
        next_edge = self.edges[next_edge_index]
        if next_edge[0] == vertex:
            return next_edge_index, next_edge[1], next_edge[0]
        else:
            return next_edge_index, next_edge[0], next_edge[1]

    def find_next_edge_c(self, edge_index, vertex):
        """ Finds next edge looping 
        couter-clockwise. 
        """
        local_index = self.vert_to_edge[vertex].index(edge_index)
        next_edge = self.vert_to_edge[vertex][local_index-1]
        return next_edge

    def is_point_in_polygon(self, point, polygon, point_2 = None):
        if point_2 == None:
            point_2 = point + np.array([1., 0.])
        intersection_count = 0
        
        for [edge, direction] in polygon:
            [b1, b2] = self.edges[edge]
            [p, ts] = line_line_intersection(self.vertices[b1], self.vertices[b2], point, point_2)
            if abs(ts[0])<1.e-6 or abs(ts[0]-1.)<1.e-6:
                return self.is_point_in_polygon(point, 
                                                polygon, 
                                                point_2 = point+np.array([random.random(), random.random()]))

            elif ts[0]+self.threshold >= 0. and ts[0]-self.threshold <=1. and ts[1]+self.threshold>=0.:
                intersection_count += 1
        
        if intersection_count%2==0:
            return False

        return True
    
    def do_remove_fractures_from_holes(self, line):
        """ Removes fracture segments if 
        they intersect a hole. 
        """                
        for hole in [self.holes[2]]:
            intersections = []
            point_2  = hole+np.array([random.random()*10., random.random()*10.])
            for (edge_index, edge) in enumerate(self.edges):
                if edge_index not in self.fracture_edges:
                    p1 = self.vertices[edge[0]]
                    p2 = self.vertices[edge[1]]
                    (intersect, ts) = line_line_intersection(p1, p2, hole, point_2)
                    if intersect != None:
                        if -1.e-6<=ts[0]<=1.+1.e-6:
                            intersections.append([abs(ts[1]), edge_index])

            [_, min_edge] = min(intersections)
            [v1, v2] = self.edges[min_edge]
            p1 = self.vertices[v1]
            p2 = self.vertices[v2]
            vec_1 = np.array(list(p2-p1)+[0.])
            vec_2 = np.array(list(hole-p1)+[0.])

            current_edge = min_edge
            if np.cross(vec_1, vec_2)[2]< 0.:
                current_vertex = v1
            else:
                current_vertex = v2

            [next_edge, next_vertex, other_vertex] = self.find_next_edge_vertex_cc(current_edge, current_vertex)
            to_be_removed = []
            while next_edge != min_edge:
                if next_edge in self.fracture_edges:
                    to_be_removed.append(next_edge)
                    [next_edge, next_vertex, other_vertex] = self.find_next_edge_vertex_cc(next_edge, other_vertex)
                else:
                    [next_edge, next_vertex, other_vertex] = self.find_next_edge_vertex_cc(next_edge, next_vertex)

            to_be_removed = list(set(to_be_removed))
            to_be_removed.sort()
            to_be_removed.reverse()
            for edge_index in to_be_removed:
                self.remove_edge(edge_index)
    
    def do_polygons(self, line):
        """ Compute all polygons in the graph.
        """
        show_poly = False
        line_split = line.split()
        if len(line_split)> 0:
            if line_split[0][:1] == "s":
                show_poly = True

        done_edges = set()
        polygons = []
        boundaries = []
        paths = []
        index_paths = []
        def direction_int(v1, v2):
            point1 = self.vertices[v1]
            point2 = self.vertices[v2]
            return (point2[0]-point1[0])*(point2[1]+point1[1])

        counter_c_edges = range(len(self.edges))
        c_edges = range(len(self.edges))
        
        while len(counter_c_edges) > 0:
            edge_index = counter_c_edges.pop(0)
            direction_sum = 0.
            current_polygon = [edge_index]
            [current_v1, current_v2] = self.edges[edge_index]
            original_v1 = current_v1
            current_path = [self.vertices[current_v1], self.vertices[current_v2]]
            current_index_path = [[edge_index, 1]]
            direction_sum += direction_int(current_v1, current_v2)

            ## Loop counter-clockwise
            next_edge = self.find_next_edge_cc(edge_index, current_v2)

            current_polygon.append(next_edge)
            [next_v1, next_v2] = self.edges[next_edge]

            if next_v1 == current_v2:
                next_vertex = next_v2
                direction_sum += direction_int(next_v1, next_v2)
                counter_c_edges.remove(next_edge)
                current_index_path.append([next_edge, 1])

            elif next_v2 == current_v2:
                next_vertex = next_v1
                direction_sum += direction_int(next_v2, next_v1)
                c_edges.remove(next_edge)
                current_index_path.append([next_edge, -1])

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
                    counter_c_edges.remove(next_edge)
                    current_index_path.append([next_edge, 1])

                elif next_v2 == next_vertex:
                    next_vertex = next_v1
                    direction_sum += direction_int(next_v2, next_v1)
                    c_edges.remove(next_edge)
                    current_index_path.append([next_edge, -1])

                else:
                    raise Exception("Problem traversing graph")

            if direction_sum < 0:
                boundaries += current_polygon
            else:
                polygons.append(current_polygon)
                paths.append(current_path)
                index_paths.append(current_index_path)

        while len(c_edges) > 0:
            edge_index = c_edges.pop(0)
            
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
                c_edges.remove(next_edge)
                current_index_path.append((next_edge, -1))

            elif next_v2 == current_v2:
                next_vertex = next_v1
                direction_sum += direction_int(next_v2, next_v1)
                counter_c_edges.remove(next_edge)
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
                    c_edges.remove(next_edge)
                    current_index_path.append((next_edge, -1))
                elif next_v2 == next_vertex:
                    next_vertex = next_v1                        
                    direction_sum += direction_int(next_v2, next_v1)
                    counter_c_edges.remove(next_edge)
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

        add_poly = True
        self.polygons = []
        self.boundaries = boundaries

        self.internal_boundaries = []

        for polygon in index_paths:
            for hole in self.holes:
                if self.is_point_in_polygon(hole, polygon):
                    add_poly = False
            if add_poly:
                self.polygons.append(polygon)
                add_poly = True
            else:
                add_poly = True
                self.internal_boundaries.append([])
                for [edge, direction] in polygon:
                    self.internal_boundaries[-1].append(edge)
                    
        for boundary_index in boundaries:
            if boundary_index in self.fracture_edges:
                self.fracture_edges.remove(boundary_index)

        if show_poly:
            fig=pylab.figure()
            ax=fig.add_subplot(111)
            patches = map(lambda x:Polygon(x, True), map(lambda x:np.array(x), paths))
            p = PatchCollection(patches,cmap=matplotlib.cm.jet,  alpha=1.0)
            
            p.set_array(np.random.rand(len(paths)))
            p.set_array(np.ones(len(paths)))
            ax.add_collection(p)
            p.set_clim([0., 1.])

            ax.autoscale_view(True, True, True)
            
            pylab.show()
    
    def do_set_fracture(self, line):
        """ Sets edge as fracture edge. 
        """
        edge_index = int(line.split()[0])
        if edge_index < 0:
            edge_index = edge_index + len(self.edges)
        self.fracture_edges.append(edge_index)
        
    def do_set_depth(self, line):
        """ Sets the reservoir depth when extruding 
        into 3D for mimpy mesh output. 
        """
        self.reservoir_depth = float(line)

    def do_mimpy_mesh(self, line):
        """ Produces a mimpy mesh from the 
        graph and stores in the filename specified.  
        """
        res_mesh = mesh.Mesh()
        res_mesh.use_face_shifted_centroid()

        file_name = line
        res_mesh.dim = 3

        for point in self.vertices:
            res_mesh.add_point(np.array([point[0], point[1], 0.]))
        
        for (index, point) in enumerate(self.vertices):
            new_index = res_mesh.add_point(np.array([point[0],point[1], self.reservoir_depth]))

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
                    res_mesh.set_face_shifted_centroid(new_face_index, centroid)

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
            res_mesh.set_face_shifted_centroid(top_face_index, centroid)

            bot_face_index = res_mesh.add_face(bot_face)
            
            normal = res_mesh.find_face_normal(bot_face_index)
            res_mesh.set_face_normal(bot_face_index, normal)
            (area, centroid) = res_mesh.find_face_centroid(bot_face_index)
            res_mesh.set_face_area(bot_face_index, area)
            res_mesh.set_face_real_centroid(bot_face_index, centroid)
            res_mesh.set_face_shifted_centroid(top_face_index, centroid)

            res_mesh.add_boundary_face(0,  bot_face_index, -1)
            res_mesh.add_boundary_face(1,  top_face_index, 1)

            current_cell.append(top_face_index)
            current_cell_orientations.append(-1)
            current_cell.append(bot_face_index)
            current_cell_orientations.append(1)
            
            new_cell_index = res_mesh.add_cell(current_cell, current_cell_orientations)
            
            res_mesh.set_cell_k(new_cell_index, np.eye(3))
            
            (volume, centroid) = res_mesh.find_volume_centroid(new_cell_index)
            res_mesh.set_cell_real_centroid(new_cell_index, centroid)
            res_mesh.set_cell_volume(new_cell_index, volume)

        for edge_index in self.boundaries:
            face_index = edge_to_face_map[edge_index]
            res_mesh.add_boundary_face(2, face_index, 1)
        
        for boundary_index in range(len(self.internal_boundaries)):
            res_mesh.add_boundary_marker(boundary_index+3, "hole"+str(boundary_index))
            for edge_index in self.internal_boundaries[boundary_index]:
                face_index = edge_to_face_map[edge_index]
                res_mesh.add_boundary_face(boundary_index+3, face_index, 1)

        

        res_mesh.build_frac_from_faces([edge_to_face_map[edge_index] for edge_index in self.fracture_edges])
            
        res_mesh.output_vtk_mesh(file_name, [res_mesh.get_cell_domain_all(),], ["DOMAIN"])
        
        #res_mesh.save_mesh(open(file_name+"saved", 'w'))
        
        pickle_file = open(file_name, 'w')
        pickle.dump(res_mesh, pickle_file)
        pickle_file.close()

    def do_EOF(self, line):
        return True

    def do_move_vertex(self, line):
        """ Move vertex to new location. 
        """
        line_split = line.split()
        index  = int(line_split[0])
        x = float(line_split[1])
        y = float(line_split[2])
    
        self.vertices[index] = np.array([x, y])

    def do_ectov(self, line):
        """ Returns the edges connected to 
        vertex number. 
        """
        vertex_index = int(line.split()[0])
        print self.vert_to_edge[vertex_index]

    def do_vctoe(self, line):
        """ Returns the vertices connected to 
        edge number. 
        """
        edge_index = int(line.split()[0])
        print self.edges[edge_index]
         

    def do_highlight_e(self, line):
        """ Highlighss edge.
        """
        edge_index = int(line.split()[0])
        self.highlight_edges=edge_index

    def add_vertex(self, point, point_threshold = None, check_threshold = True):
        """ Adds new vertex, returns 
        vertex index. 
        """
        if point_threshold:
            pass
        else:
            point_threshold = self.point_threshold
        vertex_index = -1
        if check_threshold:
            for (point_index, current_point) in enumerate(self.vertices):
                if np.linalg.norm(point-current_point)<point_threshold:
                    vertex_index = point_index
        if vertex_index < 0:
            self.vertices.append(point)
            vertex_index = len(self.vertices)-1            
            self.bounding_box[0][0] = min(self.bounding_box[0][0], point[0])
            self.bounding_box[1][0] = max(self.bounding_box[0][0], point[0])
            
            self.bounding_box[0][1] = min(self.bounding_box[0][1], point[1])
            self.bounding_box[1][1] = max(self.bounding_box[0][1], point[1])
        
        self.vertex_cache.append(vertex_index)
        if len(self.vertex_cache) > self.vertex_cache_size:
            self.vertex_cache.pop(0)

        return vertex_index

    def do_add_vertex(self, line):
        new_point = np.array(map(float, line.split()))
        new_index = self.add_vertex(new_point)
        print "new vertex", "[", new_index, "]", new_point

    def do_add_hole(self, line):
        new_point = np.array(map(float, line.split()))

        print "adding hole at [", len(self.holes),"]", new_point[0], new_point[1] 
        self.holes.append(new_point)

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
        index = 0
        for point in points:
            v2 = p2-point
            ## Check if colinear with existing edge.
            if  np.linalg.norm(v2) < self.threshold:
                raise Exception("edge length too small")
            if abs(v1.dot(v2)/np.linalg.norm(v2)-
                   np.linalg.norm(np.array(v1)))<self.threshold:
                return index
            if np.cross(np.array([v1[0], v1[1], 0.]),
                        np.array([v2[0], v2[1], 0.]))[2] > 0.:
                signs.append(1)
                all_neg = False
            else:
                signs.append(-1)
                all_pos = False
            index += 1

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
        print self.vert_to_edge[v_index]
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
        if v1<0 :
            v1 = self.vertex_cache[v1]
        if v2<0 :
            v2 = self.vertex_cache[v2]
        self.edges.append([v1, v2])

        self.update_v_to_e(len(self.edges)-1, v1, v2)
        return len(self.edges)-1

    def do_add_edge(self, line):
        [v1, v2] = map(int, line.split())
        edge_index = self.add_edge(v1, v2)
        print "new edge","[", edge_index , "]", v1, "<->", v2

    def update_v_to_e(self, edge_index, v1, v2):
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

    def do_svtov(self, line):
        """ Snaps edge index to second index location. 
        """
        split_line = line.split()
        v1_index = int(split_line[0])
        v2_index = int(split_line[1])
        
        self.vertices[v1_index] = np.array(self.vertices[v2_index])

    def set_edge(self, edge_index, v1, v2):
        """ Modifies vertices of existing 
        edge. 
        """
        [orig_v1, orig_v2] = self.edges[edge_index]

        self.vert_to_edge[orig_v1].remove(edge_index)
        self.vert_to_edge[orig_v2].remove(edge_index)

        if len(self.vert_to_edge[orig_v1]) == 2:
            self.update_vte_for_two_edges(orig_v1)

        if len(self.vert_to_edge[orig_v2]) == 2:
            self.update_vte_for_two_edges(orig_v2)

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
        to_be_merged = []
        to_be_removed =[]
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

                if (abs(param[0])<self.threshold or abs(param[0]-1.)<self.threshold) and \
                        (abs(param[1])<self.threshold or abs(param[1]-1.)<self.threshold):
                    if abs(param[0])<self.threshold:
                        merge1 = edge[0]
                    else:
                        merge1 = edge[1]
                    if abs(param[1])<self.threshold:
                        merge2 = current_edge[0]
                    else:
                        merge2 = current_edge[1]
                    to_be_merged.append([merge1, merge2])

                elif 0.-self.threshold<=param[0]<=1.+self.threshold \
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
                            new_point = self.add_vertex(intersection, 
                                                        point_threshold=self.threshold)
                            edge_mod.append([current_index, new_point])
                            segments.append([param[0], new_point])

        #for (merge1, merge2) in to_be_merged:
        #    self.merge_two_vertices(merge1, merge2)

        for current_mod in edge_mod:
            current_index = current_mod[0]
            point_index = current_mod[1]
            if  point_index != self.edges[current_index][1]:
                new_edge = self.add_edge(point_index, self.edges[current_index][1])
                p1 = self.vertices[self.edges[current_index][0]]
                p2 = self.vertices[point_index]
                if np.linalg.norm(p1-p2)<self.threshold:
                    print "HERE", p1
                self.set_edge(current_index, self.edges[current_index][0], point_index)
                if current_index in self.fracture_edges:
                    self.fracture_edges.append(new_edge)

        segments.sort()

        done_edges = set()
        done_edges.add(edge_index)
        if len(segments) > 0:
            new_segments = [segments[0]]
            for seg_index in range(1, len(segments)):
                if abs(segments[seg_index][0]-new_segments[-1][0])<self.threshold:
                    pass
                elif segments[seg_index][1]==new_segments[-1][1]:
                    pass
                else:
                    new_segments.append(segments[seg_index])
            segments = new_segments

            segments.append([1. , self.edges[edge_index][1]])
            self.set_edge(edge_index, self.edges[edge_index][0], segments[0][1])

            if edge_index in self.fracture_edges:
                add_to_frac = True
            else:
                add_to_frac = False

            for seg_index in range(len(segments)-1):
                v1 = segments[seg_index][1]
                v2 = segments[seg_index+1][1]
                if v1 != v2:
                    new_edge_index = self.add_edge(v1, v2)
                    done_edges.add(new_edge_index)
                    if add_to_frac:
                        self.fracture_edges.append(new_edge_index)

        self.remove_list_edges(to_be_removed)

        return done_edges

    def remove_overlap(self, edge_index):
        """ 
        """
        edge = self.edges[edge_index]
        point1 = self.vertices[edge[0]]
        point2 = self.vertices[edge[1]]

        # edges in graph to be modified 
        edge_mod = []

        # vertices on the primary line
        segments = []
        to_be_merged = []
        to_be_removed =[]
        for (current_index, current_edge) in enumerate(self.edges):
            if edge_index == current_index:
                pass
            else:
                current_point1 = self.vertices[current_edge[0]]
                current_point2 = self.vertices[current_edge[1]]
                
                if is_point_on_line(current_point1, point1, point2) and\
                        is_point_on_line(current_point2, point1, point2):
                    to_be_removed.append(current_index)
                    print current_index
                elif is_point_on_line(current_point1, point1, point2) and\
                        is_point_on_ray(current_point2, point1, point2):
                    if np.linalg.norm(point1-current_point2) <\
                            np.linalg.norm(point2-current_point2):
                        self.set_edge(current_index, edge[0], current_edge[1])
                    else:
                        self.set_edge(current_index, edge[1], current_edge[1])

                elif is_point_on_line(current_point2, point1, point2) and\
                        is_point_on_ray(current_point1, point1, point2):
                    if np.linalg.norm(point1-current_point1) <\
                            np.linalg.norm(point2-current_point1):
                        self.set_edge(current_index, edge[0], current_edge[0])
                    else:
                        self.set_edge(current_index, edge[1], current_edge[0])
                
        self.remove_list_edges(to_be_removed)

    def update_edge_numbering(self, edge_index):
        """ After removing an edge, updates 
        the edge numbering everywhere that is needed. 
        """
        for vertex_index in self.vert_to_edge:
            for (index,  current_edge_index) in enumerate(self.vert_to_edge[vertex_index]):
                if current_edge_index > edge_index:
                    self.vert_to_edge[vertex_index][index] -= 1

        for (index, frac_index) in enumerate(self.fracture_edges):
            if frac_index > edge_index:
                self.fracture_edges[index] -= 1

    def merge_two_vertices(self, v1, v2):
        """ Merges two vertices into one vertex. 
        Uses the midpoint between them as the new point. 
        """
        to_be_removed = []
        ## Remove edges that have v1 and v2 as vertices 
        ## since these will become degenerate
        for edge in self.vert_to_edge[v2]:
            [current_v1, current_v2] = self.edges[edge]
            if (current_v1 == v1 and current_v2 == v2) or \
                    (current_v1 == v2 and current_v2 == v1):
                to_be_removed.append(edge)

        print "merging", v1, v2, "v2e1 ->", self.vert_to_edge[v1]

        ## Remove edges that will overlap 
        ## with existing edges. 
        for edge1 in self.vert_to_edge[v1]:
            if self.edges[edge1][0] == v1:
                edge1_v2 = self.edges[edge1][1]
            else:
                edge1_v2 = self.edges[edge1][0]
            for edge2 in self.vert_to_edge[v2]:
                if self.edges[edge2][0] == v2:
                    edge2_v2 = self.edges[edge2][1]
                else:
                    edge2_v2 = self.edges[edge2][0]
                
                if edge1_v2 == edge2_v2:
                    to_be_removed.append(edge2)
        
        self.remove_list_edges(to_be_removed)

        to_be_set = []
        for edge in self.vert_to_edge[v2]:
            [current_v1, current_v2] = self.edges[edge]
            if current_v1 == v2:
                to_be_set.append([edge, v1, current_v2])
            elif current_v2 == v2:
                to_be_set.append([edge, current_v1, v1])

        
        for edge, vert1, vert2 in to_be_set:
            self.set_edge(edge, vert1, vert2)


    def do_merge_two_vertices(self, line):
        """ Merges two vertices into one vertex. 
        Uses the midpoint between them as the new point. 
        """
        line_split = line.split()
        v1 = int(line_split[0])
        v2 = int(line_split[1])
        new_point = self.merge_two_vertices(v1, v2)
        print "merged ", v1, "<=>", v2, "new point = ", new_point


    def do_merge_fracture_vertices(self, line):
        """ Merges vertices on fracture with other 
        vertices in the graph that are too close. 
        """
        merge_thresh = float(line.split()[0])

        def merge_aux():
            for f_edge in self.fracture_edges:
                [f_v1, f_v2] = self.edges[f_edge]
                f_p1 = self.vertices[f_v1]
                for o_edge in self.vert_to_edge[f_v1]:
                    if o_edge != f_edge:
                        o_vertices = list(self.edges[o_edge])
                        o_vertices.remove(f_v1)
                        o_v = o_vertices[0]
                        o_p = self.vertices[o_v]
                        if np.linalg.norm(o_p-f_p1)<merge_thresh:
                            self.merge_two_vertices(f_v1, o_v)
                            return merge_aux()
            return
                
        merge_aux()         

    def update_vte_for_two_edges(self, vertex):
        """ Updates the vertex to edge ordering 
        in case there are two edges left in the map. 
        This avoids a special problem that happens 
        removing edges from the graph. 
        """
        [edge1, edge2] = self.vert_to_edge[vertex]
        
        [e1_v1, e1_v2] = self.edges[edge1]
        [e2_v1, e2_v2] = self.edges[edge2]

        if e1_v1 == vertex:
            e1_v = e1_v2
        else:
            e1_v = e1_v1

        if e2_v1 == vertex:
            e2_v = e2_v2
        else:
            e2_v = e2_v1
            
        vector1 = self.vertices[e1_v]-self.vertices[vertex]
        vector2 = self.vertices[e2_v]-self.vertices[vertex]
        
        if np.cross(vector1, vector2)<0.:
            self.vert_to_edge[vertex] = [edge2, edge1]
            
    def remove_edge(self, edge_index):
        """ Removes edge from graph. 
        """
        [v1, v2] = self.edges[edge_index]
        self.vert_to_edge[v1].remove(edge_index)
        self.vert_to_edge[v2].remove(edge_index)
        
        if len(self.vert_to_edge[v1]) == 2:
            self.update_vte_for_two_edges(v1)

        if len(self.vert_to_edge[v2]) == 2:
            self.update_vte_for_two_edges(v2)

        self.edges.pop(edge_index)
        if edge_index in self.fracture_edges:
            self.fracture_edges.remove(edge_index)
        self.update_edge_numbering(edge_index)

    def remove_list_edges(self, edge_list):
        """ Removes a list of edges from the graph. 
        """
        to_be_removed = list(edge_list)
        while len(to_be_removed) > 0:
            current_edge = to_be_removed.pop()
            self.remove_edge(current_edge)
            for (index, edge) in enumerate(to_be_removed):
                if edge > current_edge:
                    to_be_removed[index] -= 1

    def do_remove_edge(self, line):
        """ Removes edge from graph. 
        """
        index = int(line.split()[0])
        if index < 0 :
            index = len(self.edges)+index
        self.remove_edge(index)
        
    def do_get_vertex(self, line):
        """ Returns vertex coordinates. 
        """
        try:
            print self.vertices[int(line)]
        except:
            print "No such vertex", line

    def do_find_vertex(self, line):
        """ Finds vertex closest to point. 
        """
        line_split = line.split()
        point = np.array([float(line_split[0]), float(line_split[1])])
        min = 100000
        min_vertex = -1
        for (vertex_index, vertex) in enumerate(self.vertices):
            if np.linalg.norm(vertex-point) < min:
                min = np.linalg.norm(vertex-point)
                min_vertex  = vertex_index
        print min_vertex, min

    def do_intersect_all(self, line):
        """ Intersect all edges. 
        """
        current_index = 0
        done_edges = set() 
        while current_index < len(self.edges):
            
            if current_index not in done_edges:
                done_edges.union(self.intersect_edge(current_index))

            current_index += 1
    
    def do_intersect_fractures(self, line):
        """ Intersect all fractures. 
        """
        for edge_index in self.fracture_edges:
            self.remove_overlap(edge_index)

        current_index = 0
        done_edges = set()
        to_do_edges = list(self.fracture_edges)
        while current_index < len(self.fracture_edges):            
            if self.fracture_edges[current_index] not in done_edges:
                done_edges.union(self.intersect_edge(self.fracture_edges[current_index]))

            current_index += 1

        #for current_edge in to_do_edges:
        #    if current_edge not in done_edges:
        #        done_edges.union(self.intersect_edge(current_edge))
        

    def do_intersect_edge(self, line):
        """ For a given edge number, intersect it 
        with all other edges in the graph. 
        """
        line_split = line.split()
        edge_index = int(line_split[0])
        if edge_index < 0 :
            edge_index = len(self.edges)+edge_index
        self.intersect_edge(edge_index)

    def do_load_commands(self, line):
        file = open(line)
        for line in file:
            if len(line.split())>0 and line[0]!="#":
                self.onecmd(line)

    def do_build_mesh(self, line):
        """ Builds rectangular mesh. Input:
        xstart xend, ystart, yend 
        [number of cells in x]
        [number of cells in y]
        """
        split_line = line.split()
        
        x_start = float(split_line[0])
        x_end = float(split_line[1])
        y_start = float(split_line[2])
        y_end = float(split_line[3])
        
        nx = int(split_line[4])
        ny = int(split_line[5])

        mesh_boundary = True
        if len(split_line)> 6:
            if split_line[6] == "noboundary":
                mesh_boundary = False

        ij_to_point = {}

        dx = (x_end-x_start)/nx
        dy = (y_end-y_start)/ny

        if mesh_boundary:
            for i in range(nx+1):
                for j in range(ny+1):
                    new_point = np.array([i*dx+x_start, j*dy+y_start])
                    new_point_index = self.add_vertex(new_point, check_threshold = False)
                    ij_to_point[(i, j)] = new_point_index

        else:
            for i in range(nx+1):
                for j in range(ny+1):
                    if i == 0 and j == 0:
                        pass
                    elif i == 0 and j == ny:
                        pass
                    elif i == nx and j == 0:
                        pass
                    elif i == nx and j == ny:
                        pass
                    else:
                        new_point = np.array([i*dx+x_start, j*dy+y_start])
                        new_point_index = self.add_vertex(new_point,
                                                          check_threshold = False)
                        ij_to_point[(i, j)] = new_point_index
        

        if mesh_boundary:
            for i in range(nx+1):
                for j in range(ny+1):
                    if i < nx:
                        self.add_edge(ij_to_point[(i, j)], 
                                      ij_to_point[(i+1, j)])

                    if j < ny:
                        self.add_edge(ij_to_point[(i, j)], 
                                      ij_to_point[(i, j+1)])

        else:
            for i in range(nx+1):
                for j in range(ny+1):
                    if i < nx and 0 <j <ny:
                        self.add_edge(ij_to_point[(i, j)], 
                                      ij_to_point[(i+1, j)])

                    if j < ny and 0<i <nx:
                        self.add_edge(ij_to_point[(i, j)], 
                                      ij_to_point[(i, j+1)])

    def do_build_bdm_mesh(self, line):
        """ Builds rectangular mesh with two degrees of freedom 
        per face. Input:
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
        i_mid_j_to_point = {}
        ij_mid_to_point = {}

        dx = x_dim/nx
        dy = y_dim/ny

        for i in range(nx+1):
            for j in range(ny+1):
                new_point = np.array([i*dx, j*dy])
                new_point_index = self.add_vertex(new_point)
                ij_to_point[(i, j)] = new_point_index
                
                if i < nx:
                    new_point = np.array([i*dx+.5*dx, j*dy])
                    new_point_index = self.add_vertex(new_point)
                
                    i_mid_j_to_point[(i, j)] = new_point_index
                    
                if j < ny:

                    new_point = np.array([i*dx, j*dy+.5*dy])
                    new_point_index = self.add_vertex(new_point)
                

                    ij_mid_to_point[(i, j)] = new_point_index

        for i in range(nx+1):
            for j in range(ny+1):
                
                if i < nx:
                    self.add_edge(ij_to_point[(i, j)], 
                                  i_mid_j_to_point[(i, j)])

                    self.add_edge(i_mid_j_to_point[(i, j)], 
                                  ij_to_point[(i+1, j)])
                
                if j < ny:
                    self.add_edge(ij_to_point[(i, j)], 
                                  ij_mid_to_point[(i, j)])

                    self.add_edge(ij_mid_to_point[(i, j)], 
                                  ij_to_point[(i, j+1)])
                  

def main():

    local_mesher = Mesher()

    if len(sys.argv)>1:
        try:
            local_mesher.do_load_commands(sys.argv[1])
            
        except IOError:
            print "Cannot open ", sys.argv[1]
    else:
        local_mesher.cmdloop()
    

if __name__ == "__main__":

    main()
