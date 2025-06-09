from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import time

from Packet import *
from Formular import *
from graph import *
from dijkstra import *

import copy
import random


class Satellite:
    def __init__(self, id, orb, spo, alt, phase_list, polar, sat_list, arrived_pkt_list, drop_pkt_list):
        self.sat_list = sat_list
        self.arrived_pkt_list = arrived_pkt_list
        self.drop_pkt_list = drop_pkt_list
        self.phase_list = phase_list
        self.sat_pkt_list = []

        # Network time [ms]
        self.time = 0

        # 위성 ID and satellite parameters
        self.id = id
        self.orb = orb
        self.spo = spo
        self.sat = orb * spo
        self.x = self.id // self.spo
        self.y = self.id % self.spo
        self.lat = 0
        self.lon = 0
        self.alt = alt
        self.polar = polar  # 극지방 범위

        # 라우팅 테이블 - list index: Dst [Best hop, Best weight, Second best hop, Second best weight]
        self.routing_table = []

        # Inter satellite link (ISL)
        # adjacent satellite : [intra1, intra2, inter1, inter2]
        self.adj_sat_index_list = [-1, -1, -1, -1]
        self.adj_sat_p_d_list = [-1, -1, -1, -1]
        self.intra_ISL_list = []
        self.inter_ISL_list = []
        self.intra_ISL_p_d = []
        self.inter_ISL_p_d = []

        self.count_ISL = [0, 0, 0, 0]

        # ISL queue
        self.queue_ISL_intra_1 = []
        self.queue_ISL_intra_2 = []
        self.queue_ISL_inter_1 = []
        self.queue_ISL_inter_2 = []
        self.queue_ISL = [self.queue_ISL_intra_1, self.queue_ISL_intra_2,
                          self.queue_ISL_inter_1, self.queue_ISL_inter_2]
        self.QOR = [0, 0, 0, 0]
        self.q_d_list = [0, 0, 0, 0]

        self.receive_queue = [[], [], [], []]

        self.dijkstra_spf = 0
        self.topology_graph = None
        self.org_topology_graph = None

        self.set_lla()
        self.set_adjacent_node()
        self.get_propagation_delay()

    # Setting adjacent satellite id
    def set_adjacent_node(self):
        # horizontal adjacent node
        h_adj_1 = self.id + self.spo
        h_adj_2 = self.id - self.spo
        # YDB
        if h_adj_1 >= self.sat:
            h_adj_1 = -1
        if h_adj_2 < 0:
            h_adj_2 = -1
        self.adj_sat_index_list[2] = h_adj_1
        self.adj_sat_index_list[3] = h_adj_2
        self.inter_ISL_list.append(h_adj_1)
        self.inter_ISL_list.append(h_adj_2)

        # vertical adjacent node
        v_adj_1 = self.id + 1
        v_adj_2 = self.id - 1
        if self.id // self.spo != v_adj_1 // self.spo:
            v_adj_1 -= self.spo
        if self.id // self.spo != v_adj_2 // self.spo:
            v_adj_2 += self.spo
        self.adj_sat_index_list[0] = v_adj_1
        self.adj_sat_index_list[1] = v_adj_2
        self.intra_ISL_list.append(v_adj_1)
        self.intra_ISL_list.append(v_adj_2)

    def set_lla(self):
        phasing_intra_plane = self.phase_list[0]
        phasing_inter_plane = self.phase_list[1]
        phasing_adjacent_plane = self.phase_list[2]

        self.lat = INIT_LATITUDE + self.y * phasing_intra_plane + self.x * phasing_adjacent_plane
        self.lon = self.x * phasing_inter_plane
        if self.lat >= 360:
            self.lat -= 360
        if self.lon >= 360:
            self.lon -= 360

    def get_propagation_delay(self):
        # vertical isl length
        length_v_isl = LENGTH_INTRA_ISL(self.alt, self.spo)
        p_d_v_isl = round(length_v_isl / PARAM_C)
        self.intra_ISL_p_d.append(p_d_v_isl)
        self.intra_ISL_p_d.append(p_d_v_isl)

        # horizontal isl length
        length_h_isl_1 = 0
        length_h_isl_2 = 0
        p_d_h_isl_1 = 0
        p_d_h_isl_2 = 0
        adj_sat_1 = self.inter_ISL_list[0]
        adj_sat_2 = self.inter_ISL_list[1]

        if self.is_in_polar(self.lat):  # current node in polar
            length_h_isl_1 = -1
            length_h_isl_2 = -1
            p_d_h_isl_1 = -1
            p_d_h_isl_2 = -1
        else:  # adj_node in polar
            phasing_adjacent_plane = self.phase_list[2]

            # adj_sat_1
            RI = self.phase_list[1]
            adj_sat_lat_1 = self.lat + phasing_adjacent_plane  # adj_node_1's latitude
            adj_sat_lat_2 = self.lat - phasing_adjacent_plane  # adj_node_2's latitude

            if adj_sat_1 == -1 or self.is_in_polar(adj_sat_lat_1):
                length_h_isl_1 = -1
                p_d_h_isl_1 = -1
            else:
                length_h_isl_1 = LENGTH_INTER_ISL(self.alt, RI, self.lat, adj_sat_lat_1)
                p_d_h_isl_1 = round(length_h_isl_1 / PARAM_C)

            if adj_sat_2 == -1 or self.is_in_polar(adj_sat_lat_2):
                length_h_isl_2 = -1
                p_d_h_isl_2 = -1
            else:
                length_h_isl_2 = LENGTH_INTER_ISL(self.alt, RI, self.lat, adj_sat_lat_2)
                p_d_h_isl_2 = round(length_h_isl_2 / PARAM_C)

        self.inter_ISL_p_d.append(p_d_h_isl_1)
        self.inter_ISL_p_d.append(p_d_h_isl_2)

    # 위도(lat)을 받아 현재 극지방에 있는지 확인
    def is_in_polar(self, lat):
        lat = lat % 180
        polar_range = 90 - self.polar
        polar_min = self.polar
        polar_max = 90 + polar_range

        if polar_min <= lat <= polar_max:
            return True
        return False

    # 초기 네트워크 토폴로지 설정 - 모든 ISL의 Propagation delay를 weight로 반영한 그래프 생성
    def set_topology(self, topology):
        self.topology_graph = copy.deepcopy(topology)
        self.dijkstra_spf = DijkstraSPF(self.topology_graph, self.id)
        self.update_routing_table()

    def update_routing_table(self):
        graph = copy.deepcopy(self.topology_graph)

        # # Add ISL links if available
        # if self.intra_ISL_list[0] != -1:
        #     graph.add_edge(self.id, self.intra_ISL_list[0], self.intra_ISL_p_d[0])
        # if self.intra_ISL_list[1] != -1:
        #     graph.add_edge(self.id, self.intra_ISL_list[1], self.intra_ISL_p_d[1])
        # if self.inter_ISL_list[0] != -1:
        #     graph.add_edge(self.id, self.inter_ISL_list[0], self.inter_ISL_p_d[0])
        # if self.inter_ISL_list[1] != -1:
        #     graph.add_edge(self.id, self.inter_ISL_list[1], self.inter_ISL_p_d[1])

        # Clear previous routing table
        self.routing_table = []

        # Calculate shortest paths using Dijkstra
        spf = DijkstraSPF(graph, self.id)
        for dst in range(self.sat):
            if self.id == dst:
                # No route to itself
                self.routing_table.append([])
                continue

            # Get best path and weight
            best_path = spf.get_path(dst)
            best_weight = spf.get_distance(dst)
            best_hop = best_path[1] if len(best_path) > 1 else -1  # Check for isolated nodes)
            if NEXT_HOP_MODE == 0:
                second_best_path = spf.get_path(dst)
                second_best_weight = spf.get_distance(dst)
                second_best_hop = best_path[1] if len(second_best_path) > 1 else -1  # Check for isolated nodes
            elif NEXT_HOP_MODE == 1:
                # Temporarily remove the best hop to find the second best path
                temp_graph = copy.deepcopy(graph)
                temp_graph.add_edge(self.id, best_hop, 99999999)
                temp_spf = DijkstraSPF(temp_graph, self.id)
                second_best_path = temp_spf.get_path(dst)
                second_best_weight = temp_spf.get_distance(dst)
                second_best_hop = best_path[1] if len(second_best_path) > 1 else -1  # Check for isolated nodes
            # YDB
            # elif NEXT_HOP_MODE == 2:
            #     best_hop, best_weight = self.get_next_hop_DBPR(dst)
            #     second_best_hop, second_best_weight = best_hop, best_weight
            elif NEXT_HOP_MODE == 2:
                second_best_path = spf.get_path(dst)
                second_best_weight = spf.get_distance(dst)
                second_best_hop = best_path[1] if len(second_best_path) > 1 else -1  # Check for isolated nodes
            elif NEXT_HOP_MODE == 3:
                second_best_path = spf.get_path(dst)
                second_best_weight = spf.get_distance(dst)
                second_best_hop = best_path[1] if len(second_best_path) > 1 else -1  # Check for isolated nodes

            # Append only the best path
            self.routing_table.append([best_hop, best_weight, second_best_hop, second_best_weight])

    # YDB
    def get_next_hop_DBPR(self, pkt):
        max_weight = float('-99999999')
        next_hop = -1
        Q_self = self.get_virtual_DDD(pkt.dst)

        for i, next_node_id in enumerate(self.adj_sat_index_list):
            if next_node_id == -1 or next_node_id == self.id:
                continue
            # prev 노드 고려 여부
            # if len(pkt.path) > 1:
            #     if next_node_id == pkt.path[-2]:
            #         continue

            if next_node_id == pkt.dst:
                return next_node_id, float('99999999')

            next_node = self.sat_list[next_node_id]
            Q_next = next_node.get_virtual_DDD(pkt.dst)

            w = Q_self - Q_next

            if w > max_weight:
                max_weight = w
                next_hop = next_node_id

        return next_hop, max_weight

    # YDB
    def get_virtual_DDD(self, dst_id):
        virtual_DDD = 0
        num_of_pkt_in_q = 0  # queue에 있는 pkt 수

        for isl_queue in self.queue_ISL:
            for pkt in isl_queue:
                if pkt.dst == dst_id:
                    num_of_pkt_in_q += 1

        dst_sat = self.sat_list[dst_id]
        distance = DISTANCE_CALCULATION(self, dst_sat)

        virtual_DDD += (distance / PARAM_C) * (num_of_pkt_in_q + 1)  # queue에 있는 pkt 수 + 현재 pkt

        return virtual_DDD

    def get_next_hop_PMPF(self, pkt):
        estimated_delay_list = [float('99999999') for i in range(4)]
        queueing_delay_list = self.get_queueing_delay()

        for i in range(4):
            if self.adj_sat_index_list[i] == -1:
                continue
            else:
                estimated_delay_list[i] = self.calculate_propagation_delay(self.sat_list[self.adj_sat_index_list[i]],
                                                                           self.sat_list[pkt.dst])  # 거리
                for queueing_delay in queueing_delay_list:
                    if queueing_delay >= PMPF_TH:
                        estimated_delay_list[i] = (1 - PMPF_ALPHA) * estimated_delay_list[i] + PMPF_ALPHA * (
                                    estimated_delay_list[i] + queueing_delay)
                        break

        min_delay = min(estimated_delay_list)
        min_delay_hop = self.adj_sat_index_list[estimated_delay_list.index(min_delay)]
        return min_delay_hop, min_delay

    def get_queueing_delay(self):
        queueing_delay_list = [len(self.queue_ISL_intra_1) / (PACKET_PER_MS / 4),
                               len(self.queue_ISL_intra_2) / (PACKET_PER_MS / 4),
                               len(self.queue_ISL_inter_1) / (PACKET_PER_MS / 4),
                               len(self.queue_ISL_inter_2) / (PACKET_PER_MS / 4)]

        return queueing_delay_list

    # TG 결과에 따라 패킷 수 만큼 패킷을 생성하고 next hop ISL queue에 enqueue
    def receive_generate_packet(self, pkt_list):
        # 무작위 dst 순서로 패킷을 생성함
        random_index = random.sample(range(len(pkt_list)), k=len(pkt_list))
        for dst in random_index:
            if pkt_list[dst] != 0:
                for n_pkt in range(pkt_list[dst]):
                    pkt = Packet(self.id, dst, self.time)  # (src, dst, 생성 시간)
                    pkt.state = 1
                    self.sat_pkt_list.append(pkt)
                    # Packet.num += 1
                    # 🕒 타이머 설정 (src to dst 전체 Propagation Delay의 5배)
                    total_propagation_delay = self.calculate_total_propagation_delay(pkt)
                    pkt.timer = total_propagation_delay * TIMER_FACTOR
                    pkt.org_timer = total_propagation_delay * TIMER_FACTOR
                    self.enqueue_packet(pkt)  # ISL queue에 enqueue

    def retransmit_packet(self):
        for pkt in self.sat_pkt_list:
            if pkt.state == 2 or pkt.state == -2:
                if pkt.is_ack == False:
                    if pkt.num_retx == 0:
                        pkt.state = 10
                    else:
                        pkt.state = 0
                    if pkt.num_retx <= RETX_LIMIT:
                        # print("RETX")
                        self.sat_pkt_list.remove(pkt)
                        new_pkt = Packet(pkt.src, pkt.dst, self.time)  # (src, dst, 생성 시간)
                        new_pkt.state = 1
                        new_pkt.num_retx = pkt.num_retx + 1
                        Packet.num += 1
                        # 🕒 타이머 설정 (src to dst 전체 Propagation Delay의 5배)
                        total_propagation_delay = self.calculate_total_propagation_delay(new_pkt)
                        new_pkt.timer = total_propagation_delay * TIMER_FACTOR
                        new_pkt.org_timer = total_propagation_delay * TIMER_FACTOR

                        self.sat_pkt_list.append(new_pkt)
                        self.enqueue_packet(new_pkt)  # ISL queue에 enqueue
                else:
                    if pkt.num_retx == 0:
                        pkt.state = 10
                    else:
                        pkt.state = 0
                    if pkt.num_retx <= RETX_LIMIT:
                        # print("RETX")
                        self.sat_pkt_list.remove(pkt)
                        new_pkt = Packet(pkt.src, pkt.dst, self.time)  # (src, dst, 생성 시간)
                        new_pkt.is_ack = False
                        new_pkt.state = 1
                        new_pkt.num_retx = pkt.num_retx + 1
                        Packet.num += 1
                        # 🕒 타이머 설정 (src to dst 전체 Propagation Delay의 5배)
                        total_propagation_delay = self.calculate_total_propagation_delay(new_pkt)
                        new_pkt.timer = total_propagation_delay * TIMER_FACTOR
                        new_pkt.org_timer = total_propagation_delay * TIMER_FACTOR

                        self.sat_pkt_list.append(new_pkt)
                        self.enqueue_packet(new_pkt)  # ISL queue에 enqueue

    def enqueue_packet(self, pkt):
        if pkt.timer <= 0:
            if pkt.state == 1:
                pkt.state = 2
            elif pkt.state == -1:
                pkt.state = -2
            else:
                return
        if pkt.is_ack == True:
            pkt.retx_cur -= 1
            best_hop = pkt.path[pkt.retx_cur]

            # print("id: ", pkt.id)
            # print("state: ", pkt.state)
            # print("is_ack: ", pkt.is_ack)
            # print("path: ", pkt.path)
            # print("org_path: ", pkt.org_path)
            # print("cur: ", pkt.path[pkt.retx_cur])
            # print("src: ", pkt.src)
            # print("dst: ", pkt.dst)
            # print("total_hop: ", pkt.total_hop)
            # print("timer: ", pkt.timer)
            # print("org_timer: ", pkt.org_timer)
            # print("num_retx: ", pkt.num_retx)
            # print("self.id: ", self.id)
            # input("!@!@@!@!@")
            # pkt.path.pop(0)
            # best_hop = pkt.path.pop(0)
            # if len(pkt.path) != 0:
            next_hop_queue = -1
            if self.adj_sat_index_list[0] == best_hop:
                next_hop_queue = 0
            elif self.adj_sat_index_list[1] == best_hop:
                next_hop_queue = 1
            elif self.adj_sat_index_list[2] == best_hop:
                next_hop_queue = 2
            elif self.adj_sat_index_list[3] == best_hop:
                next_hop_queue = 3
            else:
                print("@@@@@@")
                print("shit")
                print("id: ", pkt.id)
                print("state: ", pkt.state)
                print("is_ack: ", pkt.is_ack)
                print("org_path: ", pkt.org_path)
                print("path: ", pkt.path)
                print("cur: ", pkt.path[pkt.retx_cur])
                print("src: ", pkt.src)
                print("dst: ", pkt.dst)
                print("total_hop: ", pkt.total_hop)
                print("timer: ", pkt.timer)
                print("org_timer: ", pkt.org_timer)
                print("num_retx: ", pkt.num_retx)

                print("best_hop: ", best_hop)
                print("self.id: ", self.id)
                input()

                # print(next_hop_queue)
                # print(pkt.path)
                # input()

                # ACK는 큐 안 채운다고 가정하기 때문에 이 코드 말고 아래 코드로 구현
                # # next hop queue 공간이 없으면 패킷 drop
                # if len(self.queue_ISL[next_hop_queue]) > QUEUE_SIZE:
                #     pkt.drop_time = self.time
                #     pkt.state = -1
                #     # self.drop_pkt_list.append(pkt)
                # else:
                #     pkt.set_enqueue_time(self.time)
                #     self.queue_ISL[next_hop_queue].append(pkt)
                #     self.count_ISL[next_hop_queue] += 1

                pkt.set_enqueue_time(self.time)
                self.queue_ISL[next_hop_queue].append(pkt)
                self.count_ISL[next_hop_queue] += 1

        else:
            best_hop = self.routing_table[pkt.dst][0]
            w_best_hop = self.routing_table[pkt.dst][1]
            second_best_hop = self.routing_table[pkt.dst][2]
            w_second_best_hop = self.routing_table[pkt.dst][3]

            # 첫 번째 최적 경로의 큐 인덱스 결정
            bst_queue_index = -1
            if self.adj_sat_index_list[0] == best_hop:
                bst_queue_index = 0
            elif self.adj_sat_index_list[1] == best_hop:
                bst_queue_index = 1
            elif self.adj_sat_index_list[2] == best_hop:
                bst_queue_index = 2
            elif self.adj_sat_index_list[3] == best_hop:
                bst_queue_index = 3
            else:
                print("@@@@@@")
                print("shit!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print("id: ", pkt.id)
                print("best_hop: ", best_hop)
                print("self.id: ", self.id)
                print("pkt.org_path: ", pkt.org_path)
                print("state: ", pkt.state)
                print("is_ack: ", pkt.is_ack)
                print("path: ", pkt.path)
                print("src: ", pkt.src)
                print("dst: ", pkt.dst)
                print("total_hop: ", pkt.total_hop)
                print("timer: ", pkt.timer)
                print("org_timer: ", pkt.org_timer)
                print("num_retx: ", pkt.num_retx)
                input()

            # 두 번째 최적 경로의 큐 인덱스 결정
            bst2_queue_index = -1
            if self.adj_sat_index_list[0] == second_best_hop:
                bst2_queue_index = 0
            elif self.adj_sat_index_list[1] == second_best_hop:
                bst2_queue_index = 1
            elif self.adj_sat_index_list[2] == second_best_hop:
                bst2_queue_index = 2
            elif self.adj_sat_index_list[3] == second_best_hop:
                bst2_queue_index = 3
            else:
                print("@@@@@@")
                print("shit2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print("id: ", pkt.id)
                print("best_hop: ", best_hop)
                print("self.id: ", self.id)
                print("pkt.org_path: ", pkt.org_path)
                print("state: ", pkt.state)
                print("is_ack: ", pkt.is_ack)
                print("path: ", pkt.path)
                print("src: ", pkt.src)
                print("dst: ", pkt.dst)
                print("total_hop: ", pkt.total_hop)
                print("timer: ", pkt.timer)
                print("org_timer: ", pkt.org_timer)
                print("num_retx: ", pkt.num_retx)
                input()

            # 큐 대기 지연 계산
            bst_q_d = round(len(self.queue_ISL[bst_queue_index]) / PACKET_PER_MS)
            bst2_q_d = round(len(self.queue_ISL[bst2_queue_index]) / PACKET_PER_MS)

            pkt.total_hop += 1

            next_hop_queue = 0
            # DSP는 그저 다익스트라 기반으로 next hop 결정
            if NEXT_HOP_MODE == 0:
                next_hop_queue = bst_queue_index
            # OURS는 load balancing
            elif NEXT_HOP_MODE == 1:
                if len(self.queue_ISL[bst_queue_index]) < int(QUEUE_SIZE * QOR_TH[QOR_TH_idx]):
                    next_hop_queue = bst_queue_index
                else:
                    if pkt.dst == self.adj_sat_index_list[bst_queue_index]:
                        next_hop_queue = bst_queue_index
                    else:
                        next_hop_queue = bst2_queue_index

                if len(pkt.path) >= 2:
                    if pkt.path[-2] == self.adj_sat_index_list[bst_queue_index]:
                        next_hop_queue = bst2_queue_index
                    elif pkt.path[-2] == self.adj_sat_index_list[bst2_queue_index]:
                        next_hop_queue = bst_queue_index
            # YDB
            # elif NEXT_HOP_MODE == 2:
            #     next_hop_queue = bst_queue_index
            elif NEXT_HOP_MODE == 2:
                best_hop = self.get_next_hop_DBPR(pkt)[0]
                next_hop_queue = -1
                if self.adj_sat_index_list[0] == best_hop:
                    next_hop_queue = 0
                elif self.adj_sat_index_list[1] == best_hop:
                    next_hop_queue = 1
                elif self.adj_sat_index_list[2] == best_hop:
                    next_hop_queue = 2
                elif self.adj_sat_index_list[3] == best_hop:
                    next_hop_queue = 3
                # print(pkt)
                # print(pkt.id)
                # print(pkt.src)
                # print(pkt.dst)
                # print(pkt.path)
                # print(pkt.timer)
                # print(pkt.org_timer)
                # print("==========================")
                # print(self)
                # print(self.id)
                # print(best_hop)
                # print(next_hop_queue)
                # input()
            elif NEXT_HOP_MODE == 3:
                best_hop = self.get_next_hop_PMPF(pkt)[0]
                next_hop_queue = -1
                if self.adj_sat_index_list[0] == best_hop:
                    next_hop_queue = 0
                elif self.adj_sat_index_list[1] == best_hop:
                    next_hop_queue = 1
                elif self.adj_sat_index_list[2] == best_hop:
                    next_hop_queue = 2
                elif self.adj_sat_index_list[3] == best_hop:
                    next_hop_queue = 3

                # if pkt.src == 34 and pkt.dst == 78:
                #     print(pkt.path)

                # print(best_hop)
                # print(next_hop_queue)

        # next hop queue 공간이 없으면 패킷 drop
        if len(self.queue_ISL[next_hop_queue]) > QUEUE_SIZE:
            if pkt.is_ack == False:
                pkt.drop_time = self.time
                pkt.state = 2
                self.drop_pkt_list.append(pkt)
            else:
                pkt.drop_time = self.time
                pkt.state = -2
                self.drop_pkt_list.append(pkt)

        else:
            pkt.set_enqueue_time(self.time)
            # self.queue_ISL[bst_queue_index].append(pkt)
            self.queue_ISL[next_hop_queue].append(pkt)
            self.count_ISL[next_hop_queue] += 1
        # if len(pkt.path) > 1:
        #     if LENGTH_INTER_ISL(self.alt, 0, self.sat_list[pkt.path[-1]].lat, self.sat_list[pkt.path[-2]].lat) > 0:
        #         pkt.timer -= round(LENGTH_INTER_ISL(self.alt, 0, self.sat_list[pkt.path[-1]].lat, self.sat_list[pkt.path[-2]].lat) / PARAM_C)
        #         print(round(LENGTH_INTER_ISL(self.alt, 0, self.sat_list[pkt.path[-1]].lat,
        #                                      self.sat_list[pkt.path[-2]].lat) / PARAM_C))
        #     else:
        #         pkt.timer -= round(LENGTH_INTRA_ISL(self.alt, NUM_OF_SPO) / PARAM_C)
        #         print(round(LENGTH_INTRA_ISL(self.alt, NUM_OF_SPO) / PARAM_C))

    def calculate_total_propagation_delay(self, pkt):
        # 출발지에서 목적지까지의 최적 경로
        best_path = self.dijkstra_spf.get_path_between(pkt.src, pkt.dst)
        total_distance = 0
        # 경로 상의 모든 홉 사이의 Propagation Delay 합산
        for i in range(len(best_path) - 1):
            v1 = best_path[i]
            v2 = best_path[i + 1]
            if LENGTH_INTER_ISL(self.alt, 0, self.sat_list[v1].lat, self.sat_list[v2].lat) > 0:
                total_distance += LENGTH_INTER_ISL(self.alt, 0, self.sat_list[v1].lat, self.sat_list[v2].lat)
            else:
                total_distance += LENGTH_INTRA_ISL(self.alt, NUM_OF_SPO)

        # Propagation Delay 계산 (빛의 속도 기준)
        total_propagation_delay = round(total_distance / PARAM_C)
        return total_propagation_delay

    def calculate_propagation_delay(self, src, dst):
        # 출발지에서 목적지까지의 최적 경로
        best_path = src.dijkstra_spf.get_path_between(src.id, dst.id)
        total_distance = 0
        # 경로 상의 모든 홉 사이의 Propagation Delay 합산
        for i in range(len(best_path) - 1):
            v1 = best_path[i]
            v2 = best_path[i + 1]
            if LENGTH_INTER_ISL(src.alt, 0, src.sat_list[v1].lat, src.sat_list[v2].lat) > 0:
                total_distance += LENGTH_INTER_ISL(src.alt, 0, src.sat_list[v1].lat, src.sat_list[v2].lat)
            else:
                total_distance += LENGTH_INTRA_ISL(src.alt, NUM_OF_SPO)

        # Propagation Delay 계산 (빛의 속도 기준)
        total_propagation_delay = round(total_distance / PARAM_C)
        return total_propagation_delay

    def sending_packet(self):
         for pkt_ms in range(PACKET_PER_MS):
            if len(self.queue_ISL[0]) != 0:
                pkt = self.queue_ISL[0].pop(0)

                if pkt.timer <= 0:
                    if pkt.state == 1:
                        input("1")
                        pkt.state = 2
                    elif pkt.state == -1:
                        pkt.state = -2
                        input("2")
                    # continue
                if not (pkt.state == 2 or pkt.state == -2):

                    if pkt.timer <= 0:
                        print("id: ", pkt.id)
                        print("state: ", pkt.state)
                        print("is_ack: ", pkt.is_ack)
                        print("path: ", pkt.path)
                        print("cur: ", pkt.path[pkt.retx_cur])
                        print("src: ", pkt.src)
                        print("dst: ", pkt.dst)
                        print("total_hop: ", pkt.total_hop)
                        print("timer: ", pkt.timer)
                        print("org_timer: ", pkt.org_timer)
                        print("num_retx: ", pkt.num_retx)
                        input("!!")

                    next_hop = self.intra_ISL_list[0]
                    p_d = self.intra_ISL_p_d[0]
                    q_d = self.time - pkt.get_enqueue_time()
                    pkt.set_sending_next_hop(next_hop, self.time, p_d, q_d)
                    self.sat_list[next_hop].receive_pkt(pkt, 1)
            if len(self.queue_ISL[1]) != 0:
                pkt = self.queue_ISL[1].pop(0)

                if pkt.timer <= 0:
                    if pkt.state == 1:
                        pkt.state = 2
                        input("3")
                    elif pkt.state == -1:
                        pkt.state = -2
                        input("4")
                    # continue
                if not (pkt.state == 2 or pkt.state == -2):

                    if pkt.timer <= 0:
                        print("id: ", pkt.id)
                        print("state: ", pkt.state)
                        print("is_ack: ", pkt.is_ack)
                        print("path: ", pkt.path)
                        print("cur: ", pkt.path[pkt.retx_cur])
                        print("src: ", pkt.src)
                        print("dst: ", pkt.dst)
                        print("total_hop: ", pkt.total_hop)
                        print("timer: ", pkt.timer)
                        print("org_timer: ", pkt.org_timer)
                        print("num_retx: ", pkt.num_retx)
                        input("!!")

                    next_hop = self.intra_ISL_list[1]
                    p_d = self.intra_ISL_p_d[1]
                    q_d = self.time - pkt.get_enqueue_time()
                    pkt.set_sending_next_hop(next_hop, self.time, p_d, q_d)
                    self.sat_list[next_hop].receive_pkt(pkt, 0)
            if len(self.queue_ISL[2]) != 0:
                pkt = self.queue_ISL[2].pop(0)

                if pkt.timer <= 0:
                    if pkt.state == 1:
                        pkt.state = 2
                        input("5")
                    elif pkt.state == -1:
                        pkt.state = -2
                        input("6")
                    # continue
                if not (pkt.state == 2 or pkt.state == -2):

                    if pkt.timer <= 0:
                        print("id: ", pkt.id)
                        print("state: ", pkt.state)
                        print("is_ack: ", pkt.is_ack)
                        print("path: ", pkt.path)
                        print("cur: ", pkt.path[pkt.retx_cur])
                        print("src: ", pkt.src)
                        print("dst: ", pkt.dst)
                        print("total_hop: ", pkt.total_hop)
                        print("timer: ", pkt.timer)
                        print("org_timer: ", pkt.org_timer)
                        print("num_retx: ", pkt.num_retx)
                        input("!!")

                    next_hop = self.inter_ISL_list[0]
                    p_d = self.inter_ISL_p_d[0]
                    q_d = self.time - pkt.get_enqueue_time()
                    pkt.set_sending_next_hop(next_hop, self.time, p_d, q_d)
                    self.sat_list[next_hop].receive_pkt(pkt, 3)
            if len(self.queue_ISL[3]) != 0:
                pkt = self.queue_ISL[3].pop(0)

                if pkt.timer <= 0:
                    if pkt.state == 1:
                        pkt.state = 2
                        input("7")
                    elif pkt.state == -1:
                        pkt.state = -2
                        input("8")
                    # continue
                if not (pkt.state == 2 or pkt.state == -2):

                    if pkt.timer <= 0:
                        print("id: ", pkt.id)
                        print("state: ", pkt.state)
                        print("is_ack: ", pkt.is_ack)
                        print("path: ", pkt.path)
                        print("cur: ", pkt.path[pkt.retx_cur])
                        print("src: ", pkt.src)
                        print("dst: ", pkt.dst)
                        print("total_hop: ", pkt.total_hop)
                        print("timer: ", pkt.timer)
                        print("org_timer: ", pkt.org_timer)
                        print("num_retx: ", pkt.num_retx)
                        input("!!")

                    next_hop = self.inter_ISL_list[1]
                    p_d = self.inter_ISL_p_d[1]
                    q_d = self.time - pkt.get_enqueue_time()
                    pkt.set_sending_next_hop(next_hop, self.time, p_d, q_d)
                    self.sat_list[next_hop].receive_pkt(pkt, 2)

    def receive_pkt(self, pkt, isl):
        self.receive_queue[isl].append(pkt)
        # print(self.id, " is receive pkt")

    def check_receive_pkt(self):
        for receive_q in self.receive_queue:
            while len(receive_q) != 0 and receive_q[0].is_arrival(self.time):
                pkt = receive_q.pop(0)
                if pkt.timer <= 0:
                    if pkt.state == 1:
                        pkt.state = 2
                    elif pkt.state == -1:
                        pkt.state = -2
                    # continue
                if not (pkt.state == 2 or pkt.state == -2):

                    if pkt.timer <= 0:
                        print("id: ", pkt.id)
                        print("state: ", pkt.state)
                        print("is_ack: ", pkt.is_ack)
                        print("path: ", pkt.path)
                        print("cur: ", pkt.path[pkt.retx_cur])
                        print("src: ", pkt.src)
                        print("dst: ", pkt.dst)
                        print("total_hop: ", pkt.total_hop)
                        print("timer: ", pkt.timer)
                        print("org_timer: ", pkt.org_timer)
                        print("num_retx: ", pkt.num_retx)
                        input("!!")

                    if pkt.is_ack == False:
                    # if pkt.state == 1:
                        if pkt.dst == self.id:
                            self.arrived_pkt_list.append(pkt)
                            pkt.is_ack = True
                            pkt.state = 3
                            pkt.state = -1
                            # tmp = pkt.src
                            # pkt.src = pkt.dst
                            # pkt.dst = tmp
                            pkt.retx_cur = -1
                            pkt.org_path = copy.deepcopy(pkt.path)
                            # pkt.path.reverse()
                            # #
                            # print("id: ", pkt.id)
                            # print("state: ", pkt.state)
                            # print("is_ack: ", pkt.is_ack)
                            # print("path: ", pkt.path)
                            # print("org_path: ", pkt.org_path)
                            # print("cur: ", pkt.path[pkt.retx_cur])
                            # print("src: ", pkt.src)
                            # print("dst: ", pkt.dst)
                            # print("total_hop: ", pkt.total_hop)
                            # print("timer: ", pkt.timer)
                            # print("org_timer: ", pkt.org_timer)
                            # print("num_retx: ", pkt.num_retx)
                            # print("self.id: ", self.id)
                            # input("!@!@@!@!@")
                            self.enqueue_packet(pkt)
                        else:
                            self.enqueue_packet(pkt)
                    else:
                        if pkt.src == self.id:
                            pkt.state = -3
                        else:
                            self.enqueue_packet(pkt)

    # def check_receive_pkt(self):
    #     for isl in range(4):
    #         receive_q = self.receive_queue[isl]
    #         while len(receive_q) != 0 and receive_q[0].is_arrival(self.time):
    #             pkt = receive_q.pop(0)
    #
    #             # 패킷 도착 처리 (ACK는 별도로 처리)
    #             if pkt.dst == self.id:
    #                 if not pkt.is_ack:  # 일반 패킷일 때만 카운트
    #                     self.arrived_pkt_list.append(pkt)
    #
    #                 # ACK 메시지 생성 (ACK는 드롭 카운트에서 제외)
    #                 ack_pkt = Packet(pkt.dst, pkt.src, self.time)
    #                 ack_pkt.path = pkt.path[::-1]  # 도착까지의 경로를 역순으로 설정
    #                 ack_pkt.set_enqueue_time(self.time)
    #                 ack_pkt.is_ack = True  # ACK 플래그 설정
    #
    #                 # ACK를 다음 홉 큐에 추가
    #                 next_hop = ack_pkt.path[1]
    #                 queue_index = self.get_queue_index(next_hop)
    #
    #                 if queue_index == -1 or len(self.queue_ISL[queue_index]) > QUEUE_SIZE:
    #                     ack_pkt.set_drop()
    #                     self.drop_pkt_list.append(ack_pkt)
    #                 else:
    #                     self.queue_ISL[queue_index].append(ack_pkt)
    #
    #             # 도착이 아닌 경우 다음 큐로 이동
    #             else:
    #                 self.enqueue_packet(pkt)
    #
    # # 다음 홉을 찾는 함수 (Queue 인덱스 반환)
    # def get_queue_index(self, next_hop):
    #     if next_hop == self.adj_sat_index_list[0]:
    #         return 0
    #     elif next_hop == self.adj_sat_index_list[1]:
    #         return 1
    #     elif next_hop == self.adj_sat_index_list[2]:
    #         return 2
    #     elif next_hop == self.adj_sat_index_list[3]:
    #         return 3
    #     return -1

    def time_tic(self):
        self.check_receive_pkt()
        self.sending_packet()
        self.retransmit_packet()
        self.time += 1

        for i in [0, 1, 2, 3]:
            self.QOR[i] += len(self.queue_ISL[i]) / QUEUE_SIZE
            self.q_d_list[i] = len(self.queue_ISL[i]) / PACKET_PER_MS

