

import cmd 
import matplotlib.pyplot as plt
import numpy as np

plt.ion()

class helloworld(cmd.Cmd):
    
    show_edges = True
    show_vertices = True

    show_edge_numbering = True 
    show_vertex_numbering = True

    def __init__(self):
        cmd.Cmd.__init__(self)

        self.vertices = []
        self.matrix_edges = []
        self.frac_edges = []

    def do_show(self, list):
        plt.clf()
        if self.show_vertices:
            for (index, point) in enumerate(self.vertices):
                plt.scatter(point[0], point[1])
                if self.show_vertex_numbering:
                    plt.text(point[0], point[1], str(index))

        if self.show_edges:
            for (index, edge) in enumerate(self.matrix_edges):
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

    def do_add_m_edge(self, list):
        [e1, e2] = map(int, list.split())
        self.matrix_edges.append([e1, e2])
        print "new edge","[", len(self.matrix_edges)-1 , "]", e1, "<->", e2

    def do_add_f_edge(self, list):
        [e1, e2] = map(int, list.split())
        self.frac_edges.append([e1, e2])
        print "new edge","[", len(self.frac_edges)-1 , "]", e1, "<->", e2

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

        dx = nx/x_dim
        dy = ny/y_dim

        for i in range(nx+1):
            for j in range(ny+1):
                new_point = np.array([i*dx, j*dy])
                self.vertices.append(new_point)
                new_point_index = len(self.vertices)-1
                ij_to_point[(i, j)] = new_point_index

        for i in range(nx+1):
            for j in range(ny+1):
                
                if i < nx:
                    self.matrix_edges.append([ij_to_point[(i, j)], 
                                              ij_to_point[(i+1, j)]])
                
                if j < ny:
                    self.matrix_edges.append([ij_to_point[(i, j)], 
                                              ij_to_point[(i, j+1)]])
                    
        
helloworld().cmdloop()
