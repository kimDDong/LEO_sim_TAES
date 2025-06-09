from Network import *

import csv
import os
from tqdm import tqdm

class Main:
    filename = "simul_Mode-" + str(NEXT_HOP_MODE) + "_iter-" + str(NUM_ITERATION) + '_q-' + str(QUEUE_NUM) + '_rt-' + str(RETX_LIMIT) + '_timer-' + str(TIMER_FACTOR) + ".csv"

    # 파일이 이미 존재하는 경우 삭제
    # if os.path.exists(filename):
    #     os.remove(filename)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    f = open(filename, 'a' if os.path.exists(filename) else 'w', newline='')
    main_wr = csv.writer(f)

    # If the file is newly created, write the header row
    if f.tell() == 0:
        metrics = ["ppms", "avg_e2e", "avg_pd", "avg_qd", "avg_rd", "STATE_1", "STATE_2", "STATE_3", "STATE_-1",
                   "STATE_-2", "STATE_-3", "STATE_0", "STATE_10", "STATE_11"] + [f"RETX_{i}" for i in
                                                                                 range(RETX_LIMIT + 1)] + ["n_packet",
                                                                                                           "avg_hop",
                                                                                                           "avg_remain_timer",
                                                                                                           "avg_org_timer",
                                                                                                           "avg_num_retx",
                                                                                                           "org_n_packet",
                                                                                                           "dropped_packet",
                                                                                                           "drop_rate"
                                                                                                           ]
        main_wr.writerow(metrics)

    for ppms in tqdm(PACKET_GENERATION_PER_MS_LIST):
        network = Network(main_wr)
        network.simulation_proceeding(ppms)

    f.close()

if __name__ == '__main__':
    Main()
