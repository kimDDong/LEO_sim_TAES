from Network import *

import csv
import os
from tqdm import tqdm

class Main:
    filename = "simul_Mode-" + str(NEXT_HOP_MODE) + "_iter-" + str(NUM_ITERATION) + '_' + str(QOR_TH[QOR_TH_idx]) + ".csv"

    # 파일이 이미 존재하는 경우 삭제
    # if os.path.exists(filename):
    #     os.remove(filename)

    # Open the file in append mode ('a') if it exists, otherwise in write mode ('w')
    f = open(filename, 'a' if os.path.exists(filename) else 'w', newline='')
    main_wr = csv.writer(f)

    # If the file is newly created, write the header row
    if f.tell() == 0:
        main_wr.writerow(["ppms", "avg_e2e", "avg_pd", "avg_qd", "avg_rd", "arrived_pkt", "dropped_pkt", "drop_rate", "n_packet"])

    for ppms in tqdm(PACKET_GENERATION_PER_MS_LIST):
        network = Network(main_wr)
        network.simulation_proceeding(ppms)

    f.close()

if __name__ == '__main__':
    Main()
