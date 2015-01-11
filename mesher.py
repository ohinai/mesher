
import cmd 
import matplotlib.pyplot as plt
import numpy as np
import copy 

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

    def do_greet(self, line):
        print "Hello World" 

    def do_EOF(self, line):
        return True

    def do_add_vertex(self, line):
        new_point = np.array(map(float, line.split()))
        self.vertices.append(new_point)
        print "new vertex", "[", len(self.vertices)-1, "]", new_point

    def do_add_edge(self, line):
        [e1, e2] = map(int, line.split())
        self.edges.append([e1, e2])
        print "new edge","[", len(self.edges)-1 , "]", e1, "<->", e2


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
                            self.vertices.append(intersection)
                            new_point = len(self.vertices)-1
                            edge_mod.append([current_index, new_point])
                            segments.append([param[0], new_point])

        for current_mod in edge_mod:
            current_index = current_mod[0]
            point_index = current_mod[1]
            self.edges.append([self.edges[current_index][1], point_index])
            self.edges[current_index][1] = point_index
            
        segments.sort()

        segments.append([1. , self.edges[edge_index][1]])
        self.edges[edge_index][1]=segments[0][1]
        done_edges = set()
        done_edges.add(edge_index)

        for seg_index in range(len(segments)-1):
            v1  = segments[seg_index][1]
            v2 = segments[seg_index+1][1]
            self.edges.append([v1, v2])
            done_edges.add(len(self.edges)-1)
        return done_edges

    def do_remove_edge(self, line):
        """ Removes edge from graph. 
        """
        index = int(line.split()[0])
        self.edges.pop(index)
        
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
                self.vertices.append(new_point)
                new_point_index = len(self.vertices)-1
                ij_to_point[(i, j)] = new_point_index

        for i in range(nx+1):
            for j in range(ny+1):
                
                if i < nx:
                    self.edges.append([ij_to_point[(i, j)], 
                                       ij_to_point[(i+1, j)]])
                
                if j < ny:
                    self.edges.append([ij_to_point[(i, j)], 
                                       ij_to_point[(i, j+1)]])
                    
        
mesher().cmdloop()
