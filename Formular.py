import copy
import math
import random

from Params import *

# M: The number of satellite per orbit
# R: Orbital radius
def LENGTH_INTRA_ISL(altitude, M):
    return math.sqrt(2) * (PARAM_R + altitude) * math.sqrt(1 - math.cos(math.radians(360/M)))

def LENGTH_INTER_ISL(altitude, RI, lat_A, lat_B):
    R = PARAM_R + altitude
    phi = math.radians(math.fabs(lat_A - lat_B))

    return math.sqrt(2) * R * \
        math.sqrt(1 - math.cos(phi)
                  + (1 - math.cos(math.radians(RI))) * math.cos(math.radians(lat_A)) * math.cos(math.radians(lat_B)))

# Satellite id와 생성 패킷 수를 입력하면 각 dst별 패킷 수 리스트가 출력됨
def TRAFFIC_GENERATION(sat_id, pkt_num):
    dst_lst = [0 for i in range(TD_SIZE)]

    random_weights = copy.deepcopy(T_D)
    random_weights[sat_id] = 0

    for n in range(pkt_num):
        dst_i = random.choices(range(TD_SIZE), weights=random_weights)[0]
        dst_lst[dst_i] += 1

    return dst_lst

#YDB
def DISTANCE_CALCULATION(src, dst):
    # if src.id == -1:
    #     print("@@@@@@@@@@@@@@@@@@@")
    #     return 99999999999999
    lat_src = math.radians(src.lat)
    lon_src = math.radians(src.lon)
    lat_dst = math.radians(dst.lat)
    lon_dst = math.radians(dst.lon)

    dlat = lat_dst - lat_src
    dlon = lon_dst - lon_src
    
    a = math.sin(dlat / 2) ** 2 + math.cos(lat_src) * math.cos(lat_dst) * math.sin(dlon / 2) ** 2
    c = 2 * math.asin(math.sqrt(a))

    R_earth = 6371e3 # meters
    r = R_earth + src.alt

    distance = r * c
    return distance