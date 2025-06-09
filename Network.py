import Packet
import Satellite
from Params import *
from Satellite import *
from Formular import *
from Packet import *
from graph import *
from dijkstra import *
from tqdm import tqdm

import os
import numpy as np
import csv
import time
import random


class Network:
    def __init__(self, main_wr):
        self.sat_list = []
        self.arrived_pk_list = []
        self.drop_pk_list = []
        self.n_i = 0
        self.w_qd = 0
        self.ppms = 0
        self.topology_graph = Graph()
        self.main_wr = main_wr
        self.result_arr = np.zeros((NUM_OF_SAT, NUM_OF_SAT, 6))

    def set_constellation(self, n_sat, n_orb, n_spo):
        phasing_intra_plane = 360 / n_spo
        phasing_inter_plane = 180 / n_orb
        phasing_adjacent_plane = CONSTELLATION_PARAM_F * 360 / NUM_OF_SAT
        phase_list = [phasing_intra_plane, phasing_inter_plane, phasing_adjacent_plane]

        for sat_id in range(n_sat):
            sat = Satellite(sat_id, NUM_OF_ORB, NUM_OF_SPO, SAT_HEIGHT, phase_list, POLAR_LATITUDE, self.sat_list,
                            self.arrived_pk_list, self.drop_pk_list)
            self.sat_list.append(sat)

    def set_topology_graph(self):
        for sat in self.sat_list:
            for idx, delay in enumerate(sat.intra_ISL_p_d):
                if delay != -1:
                    self.topology_graph.add_edge(sat.id, sat.intra_ISL_list[idx], delay)
            for idx, delay in enumerate(sat.inter_ISL_p_d):
                if delay != -1:
                    self.topology_graph.add_edge(sat.id, sat.inter_ISL_list[idx], delay)

        if NEXT_HOP_MODE in (1, 4):
            for _ in range(W_DIJK_ITER):
                spf = [DijkstraSPF(self.topology_graph, sat.id) for sat in self.sat_list]
                for src in self.sat_list:
                    for dst in self.sat_list:
                        if src.id == dst.id:
                            continue
                        path = spf[src.id].get_path(dst.id)
                        if len(path) < 2:
                            continue
                        weight_q = (T_D[src.id] / TD_TOTAL) * (T_D[dst.id] / (TD_TOTAL - T_D[src.id]))
                        for i in range(len(path) - 1):
                            c_weight = spf[src.id].get_edge_weight(self.topology_graph, path[i], path[i + 1])
                            new_weight = c_weight + weight_q
                            self.topology_graph.add_edge(path[i], path[i + 1], new_weight)

        for sat in self.sat_list:
            sat.set_topology(self.topology_graph)

    def packet_generation(self):
        num_src = [0] * NUM_OF_SAT
        for _ in range(self.ppms):
            src_i = random.choices(range(NUM_OF_SAT), weights=T_D)[0]
            num_src[src_i] += 1

        pkt_list = [TRAFFIC_GENERATION(i, num_src[i]) for i in range(NUM_OF_SAT)]
        for sat in self.sat_list:
            sat.receive_generate_packet(pkt_list[sat.id])

        for src in range(NUM_OF_SAT):
            for dst in range(NUM_OF_SAT):
                self.result_arr[src][dst][5] += pkt_list[src][dst]

        return pkt_list

    def time_tic(self):
        for sat in self.sat_list:
            sat.time_tic()

    def check_duplicate_id(self, target_id):
        count = sum(1 for pkt in Packet.all_packets if pkt.id == target_id)
        if count > 1:
            print("⚠️ 중복된 Packet ID가 존재합니다.")

    def get_result(self):
        retx_counts = [0] * (RETX_LIMIT + 1)
        stats = {
            'STATE_1': 0, 'STATE_2': 0, 'STATE_3': 0,
            'STATE__1': 0, 'STATE__2': 0, 'STATE__3': 0,
            'STATE_0': 0, 'STATE_10': 0, 'STATE_11': 0,
            'pd_sum': 0, 'qd_sum': 0, 'rd_sum': 0,
            'hop_sum': 0, 'remain_timer_sum': 0, 'org_timer_sum': 0, 'num_retx_sum': 0,
            'org_n_packet': 0, 'dropped_packet': 0
        }

        for pkt in Packet.all_packets:
            state = pkt.state
            if state == 1: stats['STATE_1'] += 1
            elif state == 2: stats['STATE_2'] += 1
            elif state == 3: stats['STATE_3'] += 1
            elif state == -1: stats['STATE__1'] += 1
            elif state == -2: stats['STATE__2'] += 1
            elif state == -3:
                stats['STATE__3'] += 1
                stats['pd_sum'] += sum(pkt.p_d_list)
                stats['qd_sum'] += sum(pkt.q_d_list)
                stats['rd_sum'] += pkt.num_retx * pkt.org_timer
                stats['hop_sum'] += pkt.total_hop
                stats['remain_timer_sum'] += pkt.timer
                stats['org_timer_sum'] += pkt.org_timer
                stats['org_n_packet'] += 1
            elif state == 0:
                stats['STATE_0'] += 1
                stats['org_n_packet'] += 1
            elif state == 10:
                stats['STATE_10'] += 1
                stats['dropped_packet'] += 1
                stats['org_n_packet'] += 1
            elif state == 11: stats['STATE_11'] += 1

            if pkt.num_retx <= RETX_LIMIT:
                retx_counts[pkt.num_retx] += 1
                stats['num_retx_sum'] += pkt.num_retx

        s3 = stats['STATE__3']
        avg_pd = stats['pd_sum'] / s3
        avg_qd = stats['qd_sum'] / s3
        avg_rd = stats['rd_sum'] / s3
        avg_e2e = avg_pd + avg_qd + avg_rd
        avg_hop = stats['hop_sum'] / s3
        avg_remain_timer = stats['remain_timer_sum'] / s3
        avg_org_timer = stats['org_timer_sum'] / s3
        avg_num_retx = stats['num_retx_sum'] / s3
        drop_rate = stats['dropped_packet'] / s3

        return [
            self.ppms, avg_e2e, avg_pd, avg_qd, avg_rd,
            stats['STATE_1'], stats['STATE_2'], stats['STATE_3'], stats['STATE__1'], stats['STATE__2'], stats['STATE__3'],
            stats['STATE_0'], stats['STATE_10'], stats['STATE_11'],
            *retx_counts,
            len(Packet.all_packets), avg_hop, avg_remain_timer, avg_org_timer, avg_num_retx,
            stats['org_n_packet'], stats['dropped_packet'], drop_rate
        ]

    def check_packet_states(self):
        return all(pkt.state in {-3, 0, 10} for pkt in Packet.all_packets)

    def simulation_proceeding(self, ppms):
        self.set_constellation(NUM_OF_SAT, NUM_OF_ORB, NUM_OF_SPO)
        self.set_topology_graph()
        self.ppms = ppms

        print("Simulation pkt ---", self.ppms)
        for iter_ in tqdm(range(NUM_ITERATION)):
            if iter_ % PACKET_GENERATION_CYCLE == 0:
                self.packet_generation()
            self.time_tic()

            if iter_ == NUM_ITERATION - 1:
                while not self.check_packet_states():
                    self.time_tic()

        if CSV_WRITE:
            self.main_wr.writerow(self.get_result())

        filename_org = f"Mode_{NEXT_HOP_MODE}_iter-{NUM_ITERATION}_ppms-{self.ppms}"

        if CSV_WRITE:
            filename_arrived = f'pkt/pkt-{filename_org}_{QOR_TH[QOR_TH_idx]}.csv'
            with open(filename_arrived, 'w', newline='') as f_a:
                wr_a = csv.writer(f_a)
                wr_a.writerow(["Src", "Dst", "count", "E2E", "PD", "QD", "RD", "count_dropped", "total_pkt"])
                for src in range(NUM_OF_SAT):
                    for dst in range(NUM_OF_SAT):
                        if src != dst:
                            r = self.result_arr[src][dst]
                            wr_a.writerow([src, dst, r[0], r[1], r[2], r[3], r[4], r[5]])

            filename_qor = f"QOR/OQR_Mode-{NEXT_HOP_MODE}_iter-{NUM_ITERATION}_ppms-{self.ppms}_{QOR_TH[QOR_TH_idx]}.csv"
            with open(filename_qor, 'w', newline='') as f:
                wr = csv.writer(f)
                wr.writerow(["sat.id", "ISL0", "ISL1", "ISL2", "ISL3", "count_ISL0", "count_ISL1", "count_ISL2", "count_ISL3"])
                for sat in self.sat_list:
                    qors = [sat.QOR[i] / NUM_ITERATION * 100 for i in range(4)]
                    wr.writerow([sat.id] + qors + sat.count_ISL)

        Packet.num = 0
