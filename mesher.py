
import cmd 
import matplotlib.pyplot as plt
import numpy as np

plt.ion()

def line_line_intersection(a1, a2, b1, b2):
    
    ts = np.linalg.lstsq(np.array([a2-a1, -b2+b1]).T, b1-a1)[0]
    
    a_intersect = a1 + (a2-a1)*ts[0]
    b_intersect = b1 + (b2-b1)*ts[1]
    
    if np.linalg.norm(a_intersect-b_intersect) > 1.e-5:
        return (None, [-1, -1])
    
    return (a_intersect, ts)

class helloworld(cmd.Cmd):
    
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

    def do_show(self, list):
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

    def do_greet(self, list):
        print "Hello World" 

    def do_EOF(self, list):
        return True

    def do_add_vertex(self, list):
        new_point = np.array(map(float, list.split()))
        self.vertices.append(new_point)
        print "new vertex", "[", len(self.vertices)-1, "]", new_point

    def do_add_edge(self, list):
        [e1, e2] = map(int, list.split())
        self.edges.append([e1, e2])
        print "new edge","[", len(self.edges)-1 , "]", e1, "<->", e2

    def do_intersect_edge(self, list):
        """ For a given edge number, intersect it 
        with all other edges in the graph. 
        """
        list_split = list.split()
        edge_index = int(list_split[0])
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
                if intersection != None:
                    if 0<=param[0]<=1 and 0<=param[1]<=1:
                        if abs(param[0])<self.threshold:
                            edge_mod.append([current_index, edge[0]])
                        elif abs(param[0]-1)<self.threshold:
                            edge_mod.append([current_index, edge[1]])

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
        
        for seg_index in range(len(segments)-1):
            v1  = segments[seg_index][1]
            v2 = segments[seg_index+1][1]
            self.edges.append([v1, v2])
            
    def do_load_commands(self, list):
        file = open(list)
        for line in file:
            if len(line.split())>0:
                self.onecmd(line)

    def do_build_mesh(self, list):
        """ Builds rectangular mesh. Input:
        [x dimension]
        [y dimension]
        [number of cells in x]
        [number of cells in y]
        """
        split_list = list.split()
        x_dim = float(split_list[0])
        y_dim = float(split_list[1])
        
        nx = int(split_list[2])
        ny = int(split_list[3])
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
                    
        
helloworld().cmdloop()
