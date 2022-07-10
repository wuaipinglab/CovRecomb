'''
The collection set of functions
'''

import pandas as pd
import copy
import re
import os

def normalize_position(j, geno_min, geno_max):
    V_copy = copy.deepcopy(j)
    if j != "" and j != ['']:
        for v in j:
            pos = int(v.split("_")[0])
            if pos < geno_min or pos > geno_max:
                V_copy.remove(v)
    else:
        V_copy = V_copy

    return str(V_copy)


def sort_humanly(v_list):

    def tryint(s):
        try:
            return int(s)
        except ValueError:
            return s

    def str2int(v_str):
        return [tryint(sub_str) for sub_str in re.split('([0-9]+)', v_str)]

    return sorted(v_list, key=str2int, reverse=False)


def get_keys(d, value):
    return [k for k, v in d.items() if v == value]


def merge_dict(x, y):
    for k, v in x.items():
        if k in y.keys():
            y[k] += v
        else:
            y[k] = v


def creat_dir(new_dir):
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)


def A_B_star(filename):
    if "cluster" in filename.split("_")[0]:
        A_lin = filename.split("_")[1] + "*"
        if "cluster" in filename.split("_")[2]:
            B_lin = filename.split("_")[3] + "*"
        else:
            B_lin = filename.split("_")[2]
    elif "cluster" not in filename.split("_")[0]:
        A_lin = filename.split("_")[0]
        if "cluster" in filename.split("_")[1]:
            B_lin = filename.split("_")[2] + "*"
        else:
            B_lin = filename.split("_")[1]

    return A_lin, B_lin


def del_star(A_lin):
    if "*" in A_lin:
        A_lin = A_lin.strip("*")
    else:
        A_lin = A_lin

    return A_lin


def same_cluster(epi_AB, lin_cluster, cluster_dict):
    epi_AB_cope = copy.deepcopy(epi_AB)
    for i in epi_AB_cope:
        if i in lin_cluster:
            add_lin = cluster_dict[lin_cluster[i]]
            for m in add_lin:
                epi_AB.append(m)
        else:
            epi_AB = epi_AB

    return set(epi_AB)


def get_time_epi(dct, out_date, time_epi):
    return [k for (k, v) in dct.items() if out_date <= dct[k] <= time_epi]


def detect_mode(A, B, variants_all, epi, Lineage_v):
    mutation_mode = set(variants_all[epi]) & set(Lineage_v[A] + Lineage_v[B])

    return mutation_mode


def diff_days_fm(turn_mode2, df_ABC, temp_mode, temp_ABC, n, turn_mode,
                 diff_fv):
    count_n = 0
    for m in turn_mode2:
        if len(m - temp_mode) >= int(
                diff_fv):
            count_n += 1
        else:
            continue

    if count_n == len(turn_mode2):
        df_ABC = df_ABC.append(temp_ABC.loc[n])
        turn_mode.append(temp_mode)

    elif count_n < len(turn_mode2):
        turn_mode.append(temp_mode)

    return turn_mode, df_ABC


def get_time_epi(dct, out_date, time_epi):
    return [k for (k, v) in dct.items() if out_date <= dct[k] <= time_epi]


def left_right_lin(index, df_ABC):
    mode_epi = df_ABC.loc[index, "recombination_model"]
    if mode_epi.startswith("AAAA") or mode_epi.startswith(
            "BAAAA") or mode_epi.startswith("BBAAAA") or mode_epi.startswith(
                "BBBAAAA"):
        left_lin = del_star(df_ABC.loc[index, "lineage_X"])
        right_lin = del_star(df_ABC.loc[index, "lineage_Y"])
    elif mode_epi.startswith("BBBB") or mode_epi.startswith(
            "ABBBB") or mode_epi.startswith("AABBBB") or mode_epi.startswith(
                "AAABBBB"):
        left_lin = del_star(df_ABC.loc[index, "lineage_Y"])
        right_lin = del_star(df_ABC.loc[index, "lineage_X"])

    return left_lin, right_lin


def find_AB_ances(candidate_par, left_mutation, variants_all, pango_lineage,
                  left_lin):
    A_flag = 0
    epi_ances = ""
    for cp in candidate_par:
        if A_flag == 0:
            try:
                if len(set(left_mutation) - set(variants_all[cp])
                       ) <= 2 and pango_lineage[cp] == left_lin:
                    A_flag = 1
                    epi_ances = cp
                    break
            except:
                continue
        else:
            break

    return A_flag, epi_ances


def extract_AB_bk(denovo_path, file):
    df_denovo = pd.read_csv(denovo_path + file)
    if "bk2" in file:
        small_1 = df_denovo["small_1"][0]
        big_1 = df_denovo["big_1"][0]
        small_2 = df_denovo["small_2"][0]
        big_2 = df_denovo["big_2"][0]
        bk_flag = 2
        small = big = ""

    elif "bk1" in file:
        small = df_denovo["small"][0]
        big = df_denovo["big"][0]
        bk_flag = 1
        small_1 = big_1 = small_2 = big_2 = ""

    return small_1, big_1, small_2, big_2, bk_flag, small, big


def mult_times(str, n):
    must_inA = str * n
    return must_inA


def sample_lin_time(lin_A, lin_time):
    if (lin_A in lin_time):
        timeA = lin_time[lin_A]
    else:
        timeA = ''

    return timeA


def min_pairs(dic):
    if len(dic) == 0:
        return []
    min_val = min(map(lambda v: v[1], dic.items()))

    return [item for item in dic.items() if item[1] == min_val]


def bk_count(recom_pattern):
    start = recom_pattern[0]
    change = 0
    for R in recom_pattern:
        if R != start:
            change += 1
            start = R

    return change


def obtain_pattern(AB_epi, unique_A, unique_B):
    recom_pattern = ""
    for v in AB_epi:
        if v in unique_A:
            recom_pattern = recom_pattern + "X"
        elif v in unique_B:
            recom_pattern = recom_pattern + "Y"

    return recom_pattern


def parental_lineage(M):
    select_recom = M[0].split("_")
    if len(select_recom) == 2:
        lin_A_draw = select_recom[0]
        lin_B_draw = select_recom[1]
    elif len(select_recom) == 4:
        lin_A_draw = select_recom[0] + "_" + select_recom[1]
        lin_B_draw = select_recom[2] + "_" + select_recom[3]
    elif len(select_recom) == 3:
        if "cluster" in select_recom[0]:
            lin_A_draw = select_recom[0] + "_" + select_recom[1]
            lin_B_draw = select_recom[2]
        else:
            lin_A_draw = select_recom[0]
            lin_B_draw = select_recom[1] + "_" + select_recom[2]

    return lin_A_draw, lin_B_draw


def mutation_lin_unique(same_ancient):
    clean_same_ancient = []
    for mb in same_ancient:
        if "cluster" in mb:
            clean_same_ancient.append(mb.split("_")[1])
        else:
            clean_same_ancient.append(mb)

    same_num = 0
    for saan in clean_same_ancient:
        if len(saan.split(".")) < 3:
            break
        else:
            same_01 = clean_same_ancient[0].split(".")[0:3]
            if saan.split(".")[0:3] != same_01:
                break
            elif saan.split(".")[0:3] == same_01:
                same_num += 1

    return same_num
