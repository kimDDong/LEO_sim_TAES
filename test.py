from graph import *
from dijkstra import *

sat_nodes = range(16)

graph = Graph()
graph.add_edge(0, 1, 9.333)
graph.add_edge(1, 2, 9.333)
graph.add_edge(2, 3, 9.333)
graph.add_edge(4, 5, 9.333)
graph.add_edge(5, 6, 9.333)
graph.add_edge(6, 7, 9.333)
graph.add_edge(8, 9, 9.333)
graph.add_edge(9, 10, 9.333)
graph.add_edge(10, 11, 9.333)
graph.add_edge(12, 13, 9.333)
graph.add_edge(13, 14, 9.333)
graph.add_edge(14, 15, 9.333

graph.add_edge(1, 0, 9.333)
graph.add_edge(2, 1, 9.333)
graph.add_edge(3, 2, 9.333)
graph.add_edge(5, 4, 9.333)
graph.add_edge(6, 5, 9.333)
graph.add_edge(7, 6, 9.333)
graph.add_edge(9, 8, 9.333)
graph.add_edge(10, 9, 9.333)
graph.add_edge(11, 10, 9.333)
graph.add_edge(13, 12, 9.333)
graph.add_edge(14, 13, 9.333)
graph.add_edge(15, 14, 9.333)

graph.add_edge(0, 4, 1.757)
graph.add_edge(4, 8, 1.757)
graph.add_edge(8, 12, 1.757)
graph.add_edge(4, 0, 1.757)
graph.add_edge(8, 4, 1.757)
graph.add_edge(12, 8, 1.757)
graph.add_edge(1, 5, 5.004)
graph.add_edge(5, 9, 5.004)
graph.add_edge(9, 13, 5.004)
graph.add_edge(5, 1, 5.004)
graph.add_edge(9, 5, 5.004)
graph.add_edge(13, 9, 5.004)
graph.add_edge(2, 6, 7.489)
graph.add_edge(6, 10, 7.489)
graph.add_edge(10, 14, 7.489)
graph.add_edge(6, 2, 7.489)
graph.add_edge(10, 6, 7.489)
graph.add_edge(14, 10, 7.489)
graph.add_edge(3, 7, 8.83)
graph.add_edge(7, 11, 8.83)
graph.add_edge(11, 15, 8.83)
graph.add_edge(7, 3, 8.83)
graph.add_edge(11, 7, 8.83)
graph.add_edge(15, 11, 8.83)

dijkstra_list = []
for i in range(16):
    dijkstra_list.append(DijkstraSPF(graph, i))

print("%-5s %-5s" % ("label", "distance"))
for u in sat_nodes:
    print("%-5s %8f" % (u, dijkstra_list[15].get_distance(u)))

print(dijkstra_list[0].get_path(15))
print(dijkstra_list[0].get_path(15)[1:-1])
print(dijkstra_list[1].get_path(15))

route_count = [0 for i in range(16)]
print("?")
for node in dijkstra_list:
    for dst in range(16):
        path = node.get_path(dst)
        for router in path[1:-1]:
            if router == 11:
                print(dst, " : ", node.get_path(dst))
            route_count[router] += 1

print("route count: ", route_count)
# S, T, A, B, C, D, E, F, G = nodes = list("STABCDEFG")
print("999: ", dijkstra_list[15].get_path(3))

# graph = Graph()
# graph.add_edge(S, A, 4)
# graph.add_edge(S, B, 3)
# graph.add_edge(S, D, 7)
# graph.add_edge(A, C, 1)
# graph.add_edge(B, S, 3)
# graph.add_edge(B, D, 4)
# graph.add_edge(C, E, 1)
# graph.add_edge(C, D, 3)
# graph.add_edge(D, E, 1)
# graph.add_edge(D, T, 3)
# graph.add_edge(D, F, 5)
# graph.add_edge(E, G, 2)
# graph.add_edge(G, E, 2)
# graph.add_edge(G, T, 3)
# graph.add_edge(T, F, 5)

# dijkstra = DijkstraSPF(graph, S)

# print("%-5s %-5s" % ("label", "distance"))
# for u in nodes:
#     print("%-5s %8d" % (u, dijkstra.get_distance(u)))
#
# dijkstra = DijkstraSPF(graph, A)
