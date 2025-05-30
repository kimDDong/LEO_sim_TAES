from graph import *
from dijkstra import *
from Formular import *

class Test:
    def __init__(self):
        self.polar_off = True
        self.n_orb = 8
        self.n_spo = 16
        self.n_sat = self.n_orb * self.n_spo

        self.index_sat_list = range(self.n_sat)

        self.dijkstra_list = []

        self.set_graph()
        self.set_dijkstra()
        self.get_test()


    def set_graph(self):
        self.graph = Graph()

        for cur in self.index_sat_list:
            # vertical isl
            v_adj_1 = cur + 1
            if cur//self.n_spo != v_adj_1//self.n_spo:
                v_adj_1 -= self.n_spo

            v_adj_2 = cur - 1
            if cur//self.n_spo != v_adj_2//self.n_spo:
                v_adj_2 += self.n_spo

            self.graph.add_edge(cur, v_adj_1, 9.33)
            self.graph.add_edge(cur, v_adj_2, 9.33)

            # horizontal isl
            h_adj_1 = cur + self.n_spo
            if h_adj_1 > self.n_sat:
                h_adj_1 = -1
            h_adj_2 = cur - self.n_spo
            if h_adj_2 < 0:
                h_adj_2 = -1
            # print(cur, v_adj_1, v_adj_2, h_adj_1, h_adj_2)

            w_link_11 = 9.34
            w_link_34 = 7.84
            w_link_56 = 5.32
            w_link_79 = 2.68

            r = cur % self.n_spo
            if r == 0 or r == 7 or r == 8 or r == 15:
                if self.polar_off:
                    pass
                else:
                    if h_adj_1 != -1:
                        self.graph.add_edge(cur, h_adj_1, w_link_79)
                    if h_adj_2 != -1:
                        self.graph.add_edge(cur, h_adj_2, w_link_79)

            elif r == 1 or r == 6 or r == 9 or r == 14:
                if h_adj_1 != -1:
                    self.graph.add_edge(cur, h_adj_1, w_link_56)
                if h_adj_2 != -1:
                    self.graph.add_edge(cur, h_adj_2, w_link_56)

            elif r == 2 or r == 5 or r == 10 or r == 13:
                if h_adj_1 != -1:
                    self.graph.add_edge(cur, h_adj_1, w_link_34)
                if h_adj_2 != -1:
                    self.graph.add_edge(cur, h_adj_2, w_link_34)

            elif r == 3 or r == 4 or r == 11 or r == 12:
                if h_adj_1 != -1:
                    self.graph.add_edge(cur, h_adj_1, w_link_11)
                if h_adj_2 != -1:
                    self.graph.add_edge(cur, h_adj_2, w_link_11)

    def set_dijkstra(self):
        for i in range(self.n_sat):
            self.dijkstra_list.append(DijkstraSPF(self.graph, i))

    def get_test(self):
        pkt_list = []
        for i in range(self.n_sat):
            pkt_list.append(TRAFFIC_GENERATION(i))
        print(pkt_list[0])

        path_count = [0 for i in range(self.n_sat)]
        for src in range(self.n_sat):
            for dst in range(self.n_sat):
                n_dst = pkt_list[src][dst]
                path = self.dijkstra_list[src].get_path(dst)
                for n in range(n_dst):
                    for path_node in path[1:-1]:
                        path_count[path_node] += 1

        print(path_count)

        route_count = [0 for i in range(self.n_sat)]
        print("?")
        for node in self.dijkstra_list:
            for dst in range(self.n_sat):
                path = node.get_path(dst)
                for router in path[1:-1]:
                    route_count[router] += 1
        print("route count: ", route_count)

# if __name__ == '__main__':
    # test = Test()

