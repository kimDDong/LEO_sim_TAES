# Packet.py

from Params import *

class Packet:
    num = 0
    all_packets = []  # 🆕 모든 패킷 인스턴스를 저장

    def __init__(self, src, dst, t, is_ack=False):
        if is_ack == False:
            Packet.num += 1
        self.id = Packet.num
        self.size = PACKET_SIZE
        self.src = src
        self.dst = dst

        self.init_time = t
        self.enqueue_time = 0
        self.drop_time = 0

        # self.state = 0      # alive: 0, drop: -1, retransmitted: -2, arrived: 1, arrivedACK: 2
        self.state = 11      # alive: 1, drop: 2, arrived: 3, aliveACK: -1 dropACK: -2, arrivedACK: -3, retransmitted: 0, 1st_retransmitted: 10, init: 11
        self.ttl = 0

        self.eta = 0
        self.is_ack = is_ack  # 🆕 ACK 여부 플래그 추가

        self.path = [src]
        self.org_path = []
        self.p_d_list = []
        self.q_d_list = []
        self.total_hop = 0
        self.timer = 0
        self.org_timer = 0
        self.num_retx = 0
        self.retx_cur = 0

        Packet.all_packets.append(self)  # 🆕 인스턴스 등록

    def set_sending_next_hop(self, next_hop, time, p_d, q_d):
        # if not (self.state == 1 or self.state == -1):
        #     return -1
        if self.state == 1:
            if self.timer <= 0:
                if self.state == 1:
                    self.state = 2
                    # print("q")
                elif self.state == -1:
                    self.state = -2
                    # print("w")
                return
            self.path.append(next_hop)
            self.eta = time + p_d
            self.p_d_list.append(p_d)
            self.q_d_list.append(q_d)
            self.timer -= p_d
            self.timer -= q_d
        elif self.state == -1:
            if self.timer <= 0:
                if self.state == 1:
                    self.state = 2
                    # print("q")
                elif self.state == -1:
                    self.state = -2
                    # print("w")
                return
            self.eta = time + p_d
            self.p_d_list.append(p_d)
            self.q_d_list.append(q_d)
            self.timer -= p_d
            self.timer -= q_d
        # else:
        #     if self.timer <= 0:
        #         if self.state == 1:
        #             self.state = 2
        #             print("z")
        #         elif self.state == -1:
        #             self.state = -2
        #             print("x")
        #         elif self.state == 0:
        #             self.state = 0
        #
        #
        #             print("state: ", self.state)
        #             print("id: ", self.id)
        #             print("is_ack: ", self.is_ack)
        #             print("path: ", self.path)
        #             print("org_path: ", self.org_path)
        #             print("src: ", self.src)
        #             print("dst: ", self.dst)
        #             print("total_hop: ", self.total_hop)
        #             print("timer: ", self.timer)
        #             print("org_timer: ", self.org_timer)
        #             print("num_retx: ", self.num_retx)

    def time_tic(self):
        if self.state != -1:
            self.ttl += 1
        if self.ttl > PACKET_MAX_TTL:
            self.state = -1

    def set_drop(self):
        self.state = -1

    def set_enqueue_time(self, time):
        self.enqueue_time = time

    def get_enqueue_time(self):
        return self.enqueue_time

    def is_arrival(self, time):
        if self.eta <= time:
            return True
        return False

    def get_e2e(self):
        sum_p_d = 0
        sum_q_d = 0
        for pd in self.p_d_list:
            sum_p_d += pd
        for qd in self.q_d_list:
            sum_q_d += qd
        return sum_p_d + sum_q_d, sum_p_d, sum_q_d
