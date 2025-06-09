# Simulation Parameter
NUM_EXPERIMENT = 1
CSV_WRITE = True
ONLY_HIGH_DISTRIBUTION_WRITE = True
PACKET_GENERATION_CYCLE = 1
PACKET_DROP = True

RETX_LIMIT = 2

TIMER_FACTOR = 5


# DBPR
PACKET_SIZE = 512        # [B]
QUEUE_NUM = 10
QUEUE_SIZE = PACKET_SIZE * QUEUE_NUM        # [B]
CBR = [2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4]      # [Mbps]
# CBR = [2]      # [Mbps]
PACKET_GENERATION_PER_MS_LIST = [int(a*1024*1024/(PACKET_SIZE*8)) for a in CBR]


# PMPF
# CBR = [100, 300, 500]     # [Gbps]
# PACKET_GENERATION_PER_MS_LIST = [int(a*1024*1024*1024/(PACKET_SIZE*8)) for a in CBR]
NUM_EXPERIMENT = 2 #############
NUM_ITERATION = 500   # [ms]
LINK_CAPACITY = 250      # [Mbps]
SAT_HEIGHT = 570        # [km]
NUM_OF_ORB = 8
NUM_OF_SPO = 16
NUM_OF_SAT = NUM_OF_ORB * NUM_OF_SPO
PMPF_ALPHA = 0.5

PMPF_TH = 0.5/(LINK_CAPACITY/1000)



# OURS
NUM_ITERATION = 101   # [ms]

QOR_TH = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.75, 0.85]
QOR_TH_idx = 7           # 0~10
N_HOP_FACTOR = 2



# PACKET_GENERATION_PER_MS_LIST = [128*1, 128*2, 128*3, 128*4, 128*5, 128*6, 128*7, 128*8, 128*9, 128*10]
# PACKET_GENERATION_PER_MS_LIST = [128*10]
# PACKET_GENERATION_PER_MS_LIST = [750, 813, 875, 938, 1000, 1063, 1125, 1188, 1250, 1313, 1375, 1438]

# 0: [DSP] Dijkstra shortest path (DSP)
# 1: [OURS] Dijkstra[Best, Second-best] - Propagation + (Probability dst * Weight of queueing delay)
# 2: [DBPR]
# 3: [PMPF]
# 4: [OURS2]
NEXT_HOP_MODE = 0

W_QD_ROUND = 2
W_DIJK_ITER = 10
WEIGHT_Q_D = 1

# Network Parameter
PARAM_R = 6371          # [km]
PARAM_C = 300    # [km/ms] 빛의 속도
PACKET_PER_MS = 32 # 설정 근거 필요 (원래는 33 이었음)
PACKET_MAX_TTL = 1000

ISL_OFF = 1000

# Constellation param
INIT_LATITUDE = 101.25
INCLINATION = 90
CONSTELLATION_PARAM_F = 0
POLAR_LATITUDE = 88

# TLR param
TL_GREEN = 1
TL_YELLOW = 0
TL_RED = -1

TL_DETOUR_PERCENT = 0.5

# Traffic Distribution
TD_TOTAL = 11520
TD_SIZE = 8 * 16
T_D = [60, 50, 10, 20, 40 ,10 ,10 ,10 ,10 ,10 ,60 ,10 ,40 ,280 ,540 ,90,
       10, 90, 10, 10, 20, 10, 10, 10, 10, 10, 80, 40, 60, 150, 250, 90,
       10, 210, 360, 10, 10, 10, 10, 10, 10, 10, 20, 20, 90, 100, 200, 50,
       10, 290, 660, 120, 10, 10, 10, 10, 10, 10, 10, 10, 250, 200, 170, 50,
       10, 320, 500, 80, 60, 70, 30, 10, 10, 10, 10, 120, 190, 580, 270, 40,
       10, 170, 90, 20, 160, 190, 20, 10, 10, 10, 90, 150, 200, 650, 320, 40,
       10, 10, 10, 10, 150, 30, 10, 10, 10, 10, 60, 110, 90, 430, 190, 30,
       10, 290, 240, 40, 10, 10, 10, 10, 10, 30, 40, 30, 40, 10, 40, 10]

T_D_ = [60, 50, 10, 20, 40, 10, 10, 10, 90, 540, 280, 40, 10, 60, 10, 10,
       10, 90, 10, 10, 20, 10, 10, 10, 90, 250, 150, 60, 40, 80, 10, 10,
       10, 210, 360, 10, 10, 10, 10, 10, 50, 200, 100, 90, 20, 20, 10, 10,
       10, 290, 660, 120, 10, 10, 10, 10, 50, 170, 200, 250, 10, 10, 10, 10,
       10, 320, 500, 80, 60, 70, 30, 10, 40, 270, 580, 190, 120, 10, 10, 10,
       10, 170, 90, 20, 160, 190, 20, 10, 40, 320, 650, 200, 150, 90, 10, 10,
       10, 10, 10, 10, 150, 30, 10, 10, 30, 190, 430, 90, 110, 60, 10, 10,
       10, 290, 240, 40, 10, 10, 10, 10, 10, 40, 10, 40, 30, 40, 30, 10]

T_D_2 = [[60, 50, 10, 20, 40, 10, 10, 10, 90, 540, 280, 40, 10, 60, 10, 10],
       [10, 90, 10, 10, 20, 10, 10, 10, 90, 250, 150, 60, 40, 80, 10, 10],
       [10, 210, 360, 10, 10, 10, 10, 10, 50, 200, 100, 90, 20, 20, 10, 10],
       [10, 290, 660, 120, 10, 10, 10, 10, 50, 170, 200, 250, 10, 10, 10, 10],
       [10, 320, 500, 80, 60, 70, 30, 10, 40, 270, 580, 190, 120, 10, 10, 10],
       [10, 170, 90, 20, 160, 190, 20, 10, 40, 320, 650, 200, 150, 90, 10, 10],
       [10, 10, 10, 10, 150, 30, 10, 10, 30, 190, 430, 90, 110, 60, 10, 10],
       [10, 290, 240, 40, 10, 10, 10, 10, 10, 40, 10, 40, 30, 40, 30, 10]]

DAILY_VARIATION_TRAFFIC_VOLUME = [
       0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 3.0, 3.5, 8.5, 11.0, 10.5, 9.5,
       7.6, 6.5, 7.0, 8.0, 7.0, 4.6, 2.5, 1.5, 1.2, 2.0, 1.9, 1.0
]

