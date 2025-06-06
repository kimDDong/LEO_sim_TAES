#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math

class AbstractDijkstraSPF(object):

    """ Dijkstra's shortest path algorithm, abstract class. """

    def __init__(self, G, s):
        """ Calculate shortest path from s to other nodes in G. """
        self.__dist = dist = dict()
        self.__prev = prev = dict()
        visited = set()
        queue = set()

        dist[s] = 0
        prev[s] = s
        queue.add(s)

        while queue:
            u = min(queue, key=dist.get)
            for v in self.get_adjacent_nodes(G, u):
                if v in visited:
                    continue
                alt = self.get_distance(u) + self.get_edge_weight(G, u, v)
                if alt < self.get_distance(v):
                    dist[v] = alt
                    prev[v] = u
                    queue.add(v)
            queue.remove(u)
            visited.add(u)

    @staticmethod
    def get_adjacent_nodes(G, u):
        raise NotImplementedError()

    @staticmethod
    def get_edge_weight(G, u, v):
        raise NotImplementedError()

    def get_distance(self, u):
        """ Return the length of shortest path from s to u. """
        return self.__dist.get(u, math.inf)

    def get_path(self, v):
        """ Return the shortest path to v. """
        path = [v]
        while self.__prev[v] != v:
            v = self.__prev[v]
            path.append(v)
        return path[::-1]


class DijkstraSPF(AbstractDijkstraSPF):

    @staticmethod
    def get_adjacent_nodes(G, u):
        return G.get_adjacent_nodes(u)

    @staticmethod
    def get_edge_weight(G, u, v):
        return G.get_edge_weight(u, v)

    def get_distance_between(self, src, dst):
        """
        Return the distance from src to dst using the shortest path from src.
        If src is not the root of this SPF instance, raise a warning.
        """
        if src not in self._AbstractDijkstraSPF__dist:
            raise ValueError(f"src={src} is not the root of this DijkstraSPF instance.")
        if dst not in self._AbstractDijkstraSPF__dist:
            return math.inf
        return self._AbstractDijkstraSPF__dist[dst]

    def get_path_between(self, src, dst):
        """
        Return the shortest path from src to dst.
        Only allowed if src is the root of this DijkstraSPF instance.
        """
        if src not in self._AbstractDijkstraSPF__prev:
            raise ValueError(f"src={src} is not the root of this DijkstraSPF instance.")
        if dst not in self._AbstractDijkstraSPF__prev:
            return []

        path = [dst]
        while self._AbstractDijkstraSPF__prev[dst] != dst:
            dst = self._AbstractDijkstraSPF__prev[dst]
            path.append(dst)
        return path[::-1]