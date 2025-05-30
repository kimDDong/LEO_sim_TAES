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
                                # ["count_arrived", "E2E", "PD", "QD", "count_dropped", "num generate pkt]

    def set_constellation(self, n_sat, n_orb, n_spo):
        phasing_intra_plane = 360 / n_spo
        phasing_inter_plane = 180 / n_orb  # walker-star
        phasing_adjacent_plane = CONSTELLATION_PARAM_F * 360 / NUM_OF_SAT
        phase_list = [phasing_intra_plane, phasing_inter_plane, phasing_adjacent_plane]

        for sat_id in range(n_sat):
            sat = Satellite(sat_id, NUM_OF_ORB, NUM_OF_SPO, SAT_HEIGHT, phase_list, POLAR_LATITUDE, self.sat_list,
                            self.arrived_pk_list, self.drop_pk_list)
            self.sat_list.append(sat)

    def set_topology_graph(self):
        for sat in self.sat_list:
            if sat.intra_ISL_p_d[0] != -1:
                self.topology_graph.add_edge(sat.id, sat.intra_ISL_list[0], sat.intra_ISL_p_d[0])
            if sat.intra_ISL_p_d[1] != -1:
                self.topology_graph.add_edge(sat.id, sat.intra_ISL_list[1], sat.intra_ISL_p_d[1])
            if sat.inter_ISL_p_d[0] != -1:
                self.topology_graph.add_edge(sat.id, sat.inter_ISL_list[0], sat.inter_ISL_p_d[0])
            if sat.inter_ISL_p_d[1] != -1:
                self.topology_graph.add_edge(sat.id, sat.inter_ISL_list[1], sat.inter_ISL_p_d[1])

        if NEXT_HOP_MODE == 1 or NEXT_HOP_MODE == 4:
            spf = [DijkstraSPF(self.topology_graph, sat.id) for sat in self.sat_list]
            # 이전 best hop 저장 및 수렴 관련 코드 (현재 주석 처리됨)
            # prev_best_hops = [[-1 for _ in self.sat_list] for _ in self.sat_list]
            # converged = False

            # while not converged:
            for _ in range(W_DIJK_ITER):  # W_DIJK_ITER번 반복
                # converged = True  # Assume convergence

                # 최신 그래프 기반으로 SPF 업데이트
                spf = [DijkstraSPF(self.topology_graph, sat.id) for sat in self.sat_list]

                for src in self.sat_list:
                    for dst in self.sat_list:
                        if src.id == dst.id:
                            continue

                        path = spf[src.id].get_path(dst.id)
                        if len(path) < 2:
                            continue  # 유효한 경로가 없는 경우

                        best_hop = path[1]
                        # 기존 수렴 검사 관련 코드 (주석 처리됨)
                        # if prev_best_hops[src.id][dst.id] != best_hop:
                        #     converged = False
                        #     prev_best_hops[src.id][dst.id] = best_hop

                        weight_q = (T_D[src.id] / TD_TOTAL) * (T_D[dst.id] / (TD_TOTAL - T_D[src.id]))
                        for i in range(len(path) - 1):
                            c_weight = spf[src.id].get_edge_weight(self.topology_graph, path[i], path[i + 1])
                            new_weight = c_weight + weight_q
                            self.topology_graph.add_edge(path[i], path[i + 1], new_weight)

        # if NEXT_HOP_MODE == 1:
        #     for n in range(W_DIJK_ITER):
        #         spf = []
        #         for sat in self.sat_list:
        #             spf.append(DijkstraSPF(self.topology_graph, sat.id))
        #         for src in self.sat_list:
        #             for dst in self.sat_list:
        #                 if src.id == dst.id:
        #                     continue
        #                 path = spf[src.id].get_path(dst.id)
        #                 weight_q = (T_D[src.id]/TD_TOTAL)*(T_D[dst.id]/(TD_TOTAL-T_D[src.id]))
        #                 # weight_q = round(T_D[dst.id] / (TD_TOTAL - T_D[src.id]) * WEIGHT_Q_D * (T_D[src.id] / 660), 2)
        #                 # weight_q = round(T_D[dst.id] / (TD_TOTAL - T_D[src.id]) * self.w_qd * (T_D[src.id] / 660), 2)
        #                 for i in range(len(path) - 1):
        #                     c_weight = spf[src.id].get_edge_weight(self.topology_graph, path[i], path[i + 1])
        #                     self.topology_graph.add_edge(path[i], path[i + 1], c_weight + weight_q)


        for sat in self.sat_list:
            sat.set_topology(self.topology_graph)


    def packet_generation(self):
        num_src = [0 for i in range(NUM_OF_SAT)]
        # for i in range(PACKET_GENERATION_PER_MS):
        for i in range(self.ppms):
            src_i = random.choices(range(NUM_OF_SAT), weights=T_D)[0]
            num_src[src_i] += 1

        pkt_list = []
        for i in range(NUM_OF_SAT):
            pkt_list.append(TRAFFIC_GENERATION(i, num_src[i]))

        for sat in self.sat_list:
            sat.receive_generate_packet(pkt_list[sat.id])
            sat.retransmit_packet()

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
        # Arrived packet - average e2e, queueing delay, propagation delay
        e2e_sum, p_d_sum, q_d_sum, r_d_sum= 0, 0, 0, 0
        avg_e2e, avg_pd, avg_qd, avg_rd, arrived_pkt, dropped_pkt, drop_rate, n_pakcet = 0, 0, 0, 0, 0, 0, 0, 0

        if len(self.arrived_pk_list) != 0:
            print("Get result - arrived")
            pbar = tqdm(self.arrived_pk_list)
            for pk in pbar:
                e2e, p_d, q_d = pk.get_e2e()
                e2e_sum += e2e
                p_d_sum += p_d
                q_d_sum += q_d

                self.result_arr[pk.src][pk.dst][0] += 1
                self.result_arr[pk.src][pk.dst][1] += e2e
                self.result_arr[pk.src][pk.dst][2] += p_d
                self.result_arr[pk.src][pk.dst][3] += q_d

            print("Get result - dropped")
            pbar = tqdm(self.drop_pk_list)
            for pk in pbar:
                self.result_arr[pk.src][pk.dst][4] += 1

            r_d_sum = sum(
                (pkt.num_retx * pkt.org_timer)
                    for sat in self.sat_list
                    for pkt in sat.sat_pkt_list
            )

            avg_e2e = round(e2e_sum / len(self.arrived_pk_list), 2)
            avg_pd = round(p_d_sum / len(self.arrived_pk_list), 2)
            avg_qd = round(q_d_sum / len(self.arrived_pk_list), 2)
            avg_rd = round(r_d_sum / len(self.arrived_pk_list), 2)

            avg_e2e += avg_rd
            avg_e2e = round(avg_e2e, 2)

            arrived_pkt = len(self.arrived_pk_list)
            # dropped_pkt = len(self.drop_pk_list)
            # drop_rate = round(len(self.drop_pk_list) / Packet.num * 100, 2)
            n_pakcet = self.ppms*NUM_ITERATION
            dropped_pkt = n_pakcet - arrived_pkt
            drop_rate = dropped_pkt/n_pakcet
            drop_rate = round(drop_rate*100, 2)
        print(self.ppms, avg_e2e, avg_pd, avg_qd, avg_rd, arrived_pkt, dropped_pkt, drop_rate, n_pakcet)
        print(self.ppms, avg_e2e, avg_pd, avg_qd, avg_rd, len(self.arrived_pk_list), len(self.drop_pk_list), round(len(self.drop_pk_list)/Packet.num*100, 2), len(self.arrived_pk_list) + len(self.drop_pk_list))


        return [self.ppms, avg_e2e, avg_pd, avg_qd, avg_rd, len(self.arrived_pk_list), len(self.drop_pk_list), round(len(self.drop_pk_list)/Packet.num*100, 2), len(self.arrived_pk_list) + len(self.drop_pk_list)]

    def check_packet_states(self):
        # print(len(Packet.all_packets))
        for pkt in Packet.all_packets:
            if not (pkt.state == -3 or pkt.state == 0):
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
                return 0
        return 1

    # simulation Main
    def simulation_proceeding(self, ppms):
        self.set_constellation(NUM_OF_SAT, NUM_OF_ORB, NUM_OF_SPO)

        # Network Topology Graph generation
        self.set_topology_graph()
        self.ppms = ppms

        print("Simulation pkt ---", self.ppms)
        pbar = tqdm(range(NUM_ITERATION))
        for iter_ in pbar:
            # time.sleep(0.01)
            if iter_ % PACKET_GENERATION_CYCLE == 0:
                self.packet_generation()
            self.time_tic()

            if iter_ == NUM_ITERATION-1:
                while True:
                    self.time_tic()
                    simul_state = self.check_packet_states()

                    # for pkt in Packet.all_packets:
                    #     if not (pkt.state == -3 or pkt.state == 0):
                    #         if 0 <= pkt.id <= 100:
                    #             print("id: ", pkt.id)
                    #             print("state: ", pkt.state)
                    #             print("is_ack: ", pkt.is_ack)
                    #             print("path: ", pkt.path)
                    #             print("cur: ", pkt.path[pkt.retx_cur])
                    #             print("src: ", pkt.src)
                    #             print("dst: ", pkt.dst)
                    #             print("total_hop: ", pkt.total_hop)
                    #             print("timer: ", pkt.timer)
                    #             print("org_timer: ", pkt.org_timer)
                    #             print("num_retx: ", pkt.num_retx)

                    if simul_state == 1:
                        break

        if CSV_WRITE:
            self.main_wr.writerow(self.get_result())

        filename_org = "Mode_" + str(NEXT_HOP_MODE) + '_iter-' + str(NUM_ITERATION) + "_ppms-" + str(self.ppms)

        # 패킷 src/dst 별 결과 확인
        if CSV_WRITE:
            filename_arrived = 'pkt/pkt-' + filename_org + '_' + str(QOR_TH[QOR_TH_idx]) + '.csv'

            # 파일이 이미 존재하는 경우 삭제
            if os.path.exists(filename_arrived):
                os.remove(filename_arrived)

            f_a = open(filename_arrived, 'w', newline='')
            self.wr_a = csv.writer(f_a)
            self.wr_a.writerow(["Src", "Dst", "count", "E2E", "PD", "QD", "RD", "count_dropped", "total_pkt"])

        print("Write result pkt ---")
        pbar = tqdm(range(NUM_OF_SAT))
        for src in pbar:
            for dst in range(NUM_OF_SAT):
                if src == dst:
                    continue
                else:
                    result = [src, dst, self.result_arr[src][dst][0], self.result_arr[src][dst][1],
                              self.result_arr[src][dst][2], self.result_arr[src][dst][3],
                              self.result_arr[src][dst][4], self.result_arr[src][dst][5]]
                    if CSV_WRITE:
                        self.wr_a.writerow(result)
        if CSV_WRITE:
            f_a.close()

        # 각 위성 QOR 확인
        if CSV_WRITE:
            filename = "QOR/OQR_Mode-" + str(NEXT_HOP_MODE) + "_iter-" + str(NUM_ITERATION) \
                       + '_ppms-' + str(self.ppms) + '_' + str(QOR_TH[QOR_TH_idx]) + '.csv'

            # 파일이 이미 존재하는 경우 삭제
            if os.path.exists(filename):
                os.remove(filename)

            # 새 파일 생성
            f = open(filename, 'w', newline='')
            self.wr = csv.writer(f)
            self.wr.writerow(["sat.id", "ISL0", "ISL1", "ISL2", "ISL3",
                              "count_ISL0", "count_ISL1", "count_ISL2", "count_ISL3"])


        for sat in self.sat_list:
            QOR0 = sat.QOR[0] / NUM_ITERATION * 100
            QOR1 = sat.QOR[1] / NUM_ITERATION * 100
            QOR2 = sat.QOR[2] / NUM_ITERATION * 100
            QOR3 = sat.QOR[3] / NUM_ITERATION * 100
            if CSV_WRITE:
                self.wr.writerow([sat.id, QOR0, QOR1, QOR2, QOR3, sat.count_ISL[0], sat.count_ISL[1], sat.count_ISL[2],
                                  sat.count_ISL[3]])

        pbar.close()
        if CSV_WRITE:
            f.close()

        Packet.num = 0
