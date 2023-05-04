import pandas as pd
import copy
import re
import os
import numpy as np
import csv


length = 27000
fm_threshold = 0.75
len_UXY = 4
max_bk_num = 2
cor_num = 20
dirpath = "/home/soniali/Desktop/02_recom_230203/data/2023_02_08/"
sequence_path = os.path.join(dirpath, '0_raw_data', 'nextclade.tsv')
meta_path = os.path.join(dirpath, '0_raw_data', 'metadata.tsv')
geno_min,geno_max = 266,29674
parameter_souce = "Default" # "Self-define"
mini_simi_default = 0.7838
snps_path = dirpath+"/1_filtered_data/snp_norm.txt"
fv_feature_filename = "fv_clustered_enrolled.txt"

def binary_conversion(var: int):
    assert isinstance(var, int)
    if var <= 1024:
        return f'占用 {round(var / 1024, 2)} KB内存'
    else:
        return f'占用 {round(var / (1024 ** 2), 2)} MB内存'
    
    
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
    return sorted(v_list, key=lambda x:int(re.split('([0-9]+)', x)[1]), reverse=False)


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


def del_star(A_lin):
    if "*" in A_lin:
        A_lin = A_lin.split("*")[0]
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


def detect_mode(A, B, variants_all, epi, Lineage_v):
    mutation_mode = set(variants_all[epi]) & set(Lineage_v[A] + Lineage_v[B])

    return mutation_mode


def diff_days_fm(turn_mode2, df_ABC, temp_mode, temp_ABC, n, turn_mode, diff_fv):
    count_n = 0
    for m in turn_mode2:
        if len(m - temp_mode) >= int(diff_fv):
            count_n += 1
        else:
            continue
    if count_n == len(turn_mode2):
        df_ABC = df_ABC.append(temp_ABC.loc[n])
        turn_mode.append(temp_mode)
    elif count_n < len(turn_mode2):
        turn_mode.append(temp_mode)

    return turn_mode, df_ABC


def input_meta_info(df, region, country):
    region_list = []
    country_list = []
    for i in df.index:
        epi = df.loc[i,"sample_id"]
        region_list.append(region[epi])
        country_list.append(country[epi])
    df["region"], df["country"]= region_list, country_list
    return df


def summary_epi_dist(df,dirpath_recomb):
        UA_time = df.sort_values(by="collect_date" , ascending=True)
        region_list = UA_time["region"]
        month_list = []
        time_list = UA_time["collect_date"].tolist()
        for i in time_list:
            j = i[0:7]
            month_list.append(j)

        UA_time["time_month"] = month_list
        cut_time =["2020-01","2020-02","2020-03","2020-04","2020-05","2020-06","2020-07","2020-08","2020-09","2020-10","2020-11","2020-12",\
            "2021-01","2021-02","2021-03","2021-04","2021-05","2021-06","2021-07","2021-08","2021-09","2021-10","2021-11","2021-12",\
            "2022-01","2022-02","2022-03","2022-04","2022-05","2022-06","2022-07","2022-08","2022-09","2022-10","2022-11","2022-12","2023-01","2023-02"]
        
        need_time = cut_time[cut_time.index(month_list[0]):cut_time.index(month_list[-1])+1]
        nt_r_all = {}
        for nt in need_time:
            UA_time_need = UA_time[UA_time["time_month"] == nt]
            nt_region = set(UA_time_need["region"].tolist())
            
            r_all = {}
            for nt_r in nt_region:
                nt_r_num = 0
                for j in UA_time_need["region"]:
                    if j == nt_r:
                        nt_r_num +=1
                r_all[nt_r] = nt_r_num
            nt_r_all[nt] = r_all

        df_barplot = pd.DataFrame(columns=list(set(region_list)), index=need_time)
        for t in nt_r_all:
            df_barplot.loc[t] = pd.Series(nt_r_all[t])
        df_barplot.fillna(0, inplace=True)
        df_barplot = df_barplot.astype(int)

        sum_re = 0
        for i in df_barplot:
            for j in df_barplot.index:
                sum_re +=int(df_barplot.loc[j,i])

        df_barplot.to_csv(dirpath_recomb+"2_month_region_"+str(df.shape[0])+".csv",sep=',',index=True)


def get_time_epi(dct, out_date, time_epi):
    return [k for (k, v) in dct.items() if out_date <= dct[k] <= time_epi]


def calcul_bk_2(lin_A_draw,lin_B_draw,Lineage_v,epi,epiV):
    feature_SNPA = Lineage_v[lin_A_draw]
    feature_SNPB = Lineage_v[lin_B_draw]
    A_B_shared = set(feature_SNPA) & set(feature_SNPB)
    UA_mutate = (set(feature_SNPA) & set(epiV)) - set(A_B_shared)
    UB_mutate = (set(feature_SNPB) & set(epiV)) - set(A_B_shared)
    sample_special = set(epiV) - (set(feature_SNPA) | set(feature_SNPB))

    UA_mutate_unique = []
    UB_mutate_unique = []
    sample_num = pd.DataFrame(index =[epi],columns = epiV)
    
    lin_record = ""
    for j in epiV:
        if j in A_B_shared:
            sample_num.loc[epi,j] = 3
        elif j in UA_mutate:
            sample_num.loc[epi, j] = 10
            UA_mutate_unique.append(j)
            lin_record = lin_record+"X"
        elif j in UB_mutate:
            sample_num.loc[epi, j] = -10
            UB_mutate_unique.append(j)
            lin_record = lin_record+"Y"
        elif j in sample_special:
            sample_num.loc[epi, j] = -3
    sample_num_2 = sample_num.fillna(0)
    
    return sample_num_2,lin_record,UA_mutate_unique,UB_mutate_unique


def left_right_lin(index, df_ABC):
    mode_epi = df_ABC.loc[index, "mutation_pattern"]
    if mode_epi.startswith("XXXX") or mode_epi.startswith("YXXXX") or mode_epi.startswith("YYXXXX") or mode_epi.startswith("YYYXXXX"):
        left_lin = df_ABC.loc[index, "lineage_X"]
        right_lin = df_ABC.loc[index, "lineage_Y"]
    elif mode_epi.startswith("YYYY") or mode_epi.startswith("XYYYY") or mode_epi.startswith("XXYYYY") or mode_epi.startswith("XXXYYYY"):
        left_lin = df_ABC.loc[index, "lineage_Y"]
        right_lin = df_ABC.loc[index, "lineage_X"]

    return left_lin, right_lin


def find_AB_ances(candidate_par, left_mutation, variants_all, pango_lineage, left_lin):
    A_flag = 0
    epi_ances = ""
    for cp in candidate_par:
        if A_flag == 0:
            try:
                if len(set(left_mutation) - set(variants_all[cp])) <= 2 and pango_lineage[cp] == left_lin:
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
    most_one = [item for item in dic.items() if item[1] == min_val]
    
    if len(most_one) == 1:
        return most_one
    elif len(most_one) >= 2:
        return [most_one[0]]


def bk_count(recom_pattern):
    start = recom_pattern[0]
    change = 0
    for R in recom_pattern:
        if R != start:
            change += 1
            start = R

    return change


def write_new_file(file, variants_all):
    with open(file, 'w') as f:
        for i in variants_all:
            f.write(i+':'+','.join(variants_all[i])+'\n')


def obtain_pattern(AB_epi, unique_A, unique_B):
    recom_pattern = ""
    for v in AB_epi:
        if v in unique_A:
            recom_pattern = recom_pattern + "X"
        elif v in unique_B:
            recom_pattern = recom_pattern + "Y"

    return recom_pattern


def df_score_XY(lin_X, lin_Y, df_score):
    if "*" not in lin_X:
        df_score["lineage_X"] = lin_X
    else:
        df_score["lineage_X"] = lin_X.split("*")[0]

    if "*" not in lin_Y:
        df_score["lineage_Y"] = lin_Y
    else:
        df_score["lineage_Y"] = lin_Y.split("*")[0]
    return df_score


def find_bk(df_epi_info_in, epi, variants_all, path_breakloc):
    def get_epi_info(df_epi_info_in,string):
        epi_info = df_epi_info_in.loc[df_epi_info_in.sample_id == epi, string].iloc[0]
        try: 
            epi_info2 = epi_info.split("/")
            return epi_info2
        except:
            return []
    
    lin_X = df_epi_info_in.loc[df_epi_info_in.sample_id == epi, "lineage_X"].iloc[0]
    lin_Y = df_epi_info_in.loc[df_epi_info_in.sample_id == epi, "lineage_Y"].iloc[0]

    UA_mutate = get_epi_info(df_epi_info_in,"X_mutations")
    UB_mutate = get_epi_info(df_epi_info_in,"Y_mutations")
    sample_special = get_epi_info(df_epi_info_in,"denovo_mutations")
    A_B_shared = get_epi_info(df_epi_info_in,"shared_mutations")
    feature_pattern = df_epi_info_in.loc[df_epi_info_in.sample_id == epi, "mutation_pattern"].iloc[0]
    epiV = sort_humanly(variants_all[epi])

    sample_num = pd.DataFrame(index=[epi], columns=epiV)

    for j in epiV:
        if j in A_B_shared:
            sample_num.loc[epi, j] = 3
        elif j in UA_mutate:
            sample_num.loc[epi, j] = 10
        elif j in UB_mutate:
            sample_num.loc[epi, j] = -10
        elif j in sample_special:
            sample_num.loc[epi, j] = -3

    df_score = sample_num
    variants = list(df_score)
    for i in variants:
        if (int(df_score[i]) == -10) or (int(df_score[i]) == 10):
            start = int(df_score[i])
            break

    change = 0
    for i in variants:
        if int(df_score[i]) == start:
            break_small = i.strip().split("_")[0]
        if (int(df_score[i]) != start) and (int(df_score[i]) != 3) and (int(df_score[i]) != -3):
            change += 1
            start = int(df_score[i])
            break_big = i.strip().split("_")[0]
            if change == 0:
                continue
            if change == 1:
                break_total = str(break_small) + "," + str(break_big)
            if change == 2:
                break_total = str(break_total) + "," + str(break_small) + "," + str(break_big)

    if (change == 1) and (len(break_total.split(",")) == 2):
        epi = epi.replace("/", "_")
        df_score["small"] = int(re.findall(r"\d+",break_total.split(",")[0])[0])
        df_score["big"] = int(re.findall(r"\d+",break_total.split(",")[1])[0])
        df_score = df_score_XY(lin_X, lin_Y, df_score)
        df_score["lineage_X"] = lin_X.strip().split("_" + r".*" + "free")[0].split("_")[0]
        df_score["mutation_pattern"] = feature_pattern
        df_score2 = df_score[["small", "big", "lineage_X", "lineage_Y", "mutation_pattern"]]
        df_score2.to_csv(path_breakloc + "bk" + str(change) + "_" + epi+".csv", sep=',', index=True)

    if (change == 2) and (len(break_total.split(",")) == 4):
        epi = epi.replace("/", "_")
        df_score["small_1"] = int(re.findall(r"\d+",break_total.split(",")[0])[0])
        df_score["big_1"] = int(re.findall(r"\d+",break_total.split(",")[1])[0])
        df_score["small_2"] = int(re.findall(r"\d+",break_total.split(",")[2])[0])
        df_score["big_2"] = int(re.findall(r"\d+",break_total.split(",")[3])[0])
        df_score = df_score_XY(lin_X, lin_Y, df_score)

        df_score2 = df_score[["small_1", "big_1", "small_2", "big_2", "lineage_X", "lineage_Y"]]
        df_score2.to_csv(path_breakloc + "bk" + str(change) + "_" + epi+".csv", sep=',', index=True)


def calcul_bk(lin_A_draw, lin_B_draw, Lineage_v, epiV):
    feature_SNPA = Lineage_v[lin_A_draw]
    feature_SNPB = Lineage_v[lin_B_draw]
    A_B_shared = set(feature_SNPA) & set(feature_SNPB)
    UA_mutate = (set(feature_SNPA) & set(epiV)) - set(A_B_shared)
    UB_mutate = (set(feature_SNPB) & set(epiV)) - set(A_B_shared)
    sample_special = set(epiV) - (set(feature_SNPA) | set(feature_SNPB))

    UA_mutate_unique = []
    UB_mutate_unique = []
    shared_mut = []
    denovo_mut = []

    lin_record = ""
    for j in epiV:
        if j in A_B_shared:
            shared_mut.append(j)
        elif j in UA_mutate:
            UA_mutate_unique.append(j)
            lin_record = lin_record + "X"
        elif j in UB_mutate:
            UB_mutate_unique.append(j)
            lin_record = lin_record + "Y"
        elif j in sample_special:
            denovo_mut.append(j)

    return lin_record, UA_mutate_unique, UB_mutate_unique, shared_mut, denovo_mut


def calculate_A_B_K(change,recom_pattern,AB_epi,unique_A,unique_B,Lineage_v,A,B):
    if change == 1 and recom_pattern[0] == "X":
        for v in AB_epi:
            if v in unique_A:
                A_mut = v
            elif v in unique_B:
                B_mut = v
                break
        A_K = sort_humanly(Lineage_v[A]).index(A_mut)+1
        B_K = len(Lineage_v[B]) - sort_humanly(Lineage_v[B]).index(B_mut)
    elif change == 1 and recom_pattern[0] == "Y":
        for v in AB_epi:
            if v in unique_B:
                B_mut = v
            elif v in unique_A:
                A_mut = v
                break
        B_K = sort_humanly(Lineage_v[B]).index(B_mut)+1
        A_K = len(Lineage_v[A]) - sort_humanly(Lineage_v[A]).index(A_mut)

    elif change == 2 and recom_pattern[0] == "X":
        temp_str = 0
        for v in AB_epi:
            if v in unique_A:
                temp_str+=1
                A_mut1 = v
            elif v in unique_B:
                B_mut1 = v
                temp_str+=1
                break
        for v in AB_epi[temp_str:]:
            if v in unique_B:
                B_mut2 = v
            elif v in unique_A:
                A_mut2 = v
                break
        A_K = (sort_humanly(Lineage_v[A]).index(A_mut1)+1) +(len(Lineage_v[A]) - sort_humanly(Lineage_v[A]).index(A_mut2))
        B_K = sort_humanly(Lineage_v[B]).index(B_mut2)+1-sort_humanly(Lineage_v[B]).index(B_mut1)
    elif change == 2 and recom_pattern[0] == "Y":
        temp_str = 0
        for v in AB_epi:
            if v in unique_B:
                temp_str+=1
                B_mut1 = v
            elif v in unique_A:
                A_mut1 = v
                temp_str+=1
                break
        for v in AB_epi[temp_str:]:
            if v in unique_A:
                A_mut2 = v
            elif v in unique_B:
                B_mut2 = v
                break
        B_K = (sort_humanly(Lineage_v[B]).index(B_mut1)+1) +(len(Lineage_v[B]) - sort_humanly(Lineage_v[B]).index(B_mut2))
        A_K = sort_humanly(Lineage_v[A]).index(A_mut2)+1-sort_humanly(Lineage_v[A]).index(A_mut1)
    
    return  A_K,B_K

'''
The core algorithm of the CovRecomb method to detect recombinants by assigning feature mutations.
'''

def recombination_detection(pango_lineage, len_UXY, max_bk_num, must_inA, must_inB, lin_list,\
    Strain_list_snp, collect_date, variants_all, feature_mutations, Lineage_v, mutaions_num, output_file,output_file_allp):
    
    from scipy.stats import hypergeom
    from statsmodels.sandbox.stats.multicomp import multipletests
    from function_set import bk_count, obtain_pattern,sort_humanly, bk_count,calcul_bk

    num = 0
    for epi in Strain_list_snp:
        num += 1
        if num % 1000 == 0:
            print("Processing {}:    {}/{}    proportion: {}%".format(os.getpid(), num, len(Strain_list_snp), round(num/len(Strain_list_snp),8)*100))
        try:
            epi_record = {}
            epiV = variants_all[epi]
            epi_feat = len(set(epiV) & set(feature_mutations)) # n: 样本突变个数与其LSPM矩阵中纳入候选谱系(如1000个谱系)的特征突变的交集个数；

            if epi_feat < 2 * len_UXY:
                continue
            else:
                # P-value for Non-recombination
                for lin_A in lin_list:
                    all_AA = len(Lineage_v[lin_A]) # K：正在分析的谱系对（如X+X），其特征突变嵌合在基因组上的突变个数，对于X+X，即X的特征突变个数。
                    all_AA_epi = len(set(Lineage_v[lin_A]) & set(epiV)) # k：样本突变个数与正在分析的谱系对（如X+Y）特征突变集和的交集个数
                    pVal = hypergeom.sf(all_AA_epi - 1, mutaions_num, all_AA, epi_feat)
                    epi_record[str(lin_A) + "_" + str(lin_A)] = pVal

                # the least p-value for the Non-recombinant
                min_AA = min(epi_record, key = epi_record.get)

                # P-value for Recombinant (A+B/A+B+A)
                A_already = []
                for A in lin_list:
                    A_already.append(A)
                    A_epi = set(Lineage_v[A]) & set(epiV)
                    if len(A_epi) < len_UXY:
                        continue
                    else:
                        afterA_linB = set(lin_list) - set(A_already)
                        for B in afterA_linB:
                            B_epi = set(Lineage_v[B]) & set(epiV)
                            unique_A = A_epi - B_epi
                            unique_B = B_epi - A_epi
                            if len(B_epi) < len_UXY or len(unique_A) < len_UXY or len(unique_B) < len_UXY:
                                continue
                            else:
                                union_AB_set = set(Lineage_v[A]) ^ set(Lineage_v[B])
                                AB_epi = sort_humanly(list(union_AB_set & set(epiV)))
                                recom_pattern = obtain_pattern(AB_epi, unique_A, unique_B)
                                if (must_inA not in recom_pattern) or (must_inB not in recom_pattern):
                                    continue
                                else:
                                    change = bk_count(recom_pattern)
                                    if change > max_bk_num:
                                        continue
                                    else:
                                        all_AB = len(set(Lineage_v[A]) | set(Lineage_v[B]))  # K：正在分析的谱系对（如X+Y），其特征突变嵌合在基因组上的突变个数，对于X+Y，特征突变去重后的个数。 #去重
                                        all_AB_epi = len(set(set(Lineage_v[A]) | set(Lineage_v[B])) & set(epiV)) # k：样本突变个数与正在分析的谱系对（如X+Y）特征突变集和的交集个数
                                        pVal = hypergeom.sf(all_AB_epi - 1, mutaions_num, all_AB, epi_feat)
                                        epi_record[str(A) + "_" + str(B)] = pVal

                raw_pvals = list(epi_record.values())
                rejected, p_adjusted, _, alpha_corrected = multipletests(raw_pvals, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)

                lin_adjP = {}
                for p in range(len(p_adjusted)):
                    lin_pair = list(epi_record.keys())[p]
                    lin_adjP[lin_pair] = p_adjusted[p]
                
                min_adjp_pair = min(lin_adjP, key = lin_adjP.get)
                if min_adjp_pair != min_AA and lin_adjP[min_adjp_pair] < 0.05: # 校正后p值显著，且排名top1，且其p值比PminAA小
                    epiV = sort_humanly(epiV)
                    lin_A_draw, lin_B_draw = min_adjp_pair.split("_")[0],min_adjp_pair.split("_")[1]
                    lin_record, UA_mutate_unique, UB_mutate_unique, shared_mut, denovo_mut = calcul_bk(lin_A_draw, lin_B_draw, Lineage_v, epiV)
                    with open(output_file, "a+") as file_epi:
                        file_epi.write(epi + "," + collect_date[epi] + "," + pango_lineage[epi]+","+lin_A_draw + "," + lin_B_draw + "," + lin_record + "," + "/".join(UA_mutate_unique) + "," + "/".join(UB_mutate_unique) + "," + "/".join(shared_mut) + "," + "/".join(denovo_mut) + "\n")
                else:
                    continue
                
        except:
            continue


def cluster_similar_lineages(fv_path, len_UXY, fv_feature_filename):
    Lineage_v = {}
    lineage_n = {}
    linA_list = []
    with open(fv_path, 'r') as f:
        for i in csv.reader(f): # 2302
            if (i[0].startswith("X")) or (i[0] == "Unassigned") or (i[0] == "None") or (len(i[2:]) < len_UXY) or (int(i[1]) <= 100):
                continue
            else:
                linA_list.append(i[0])
                lineage_n[i[0]] = i[1]
                Lineage_v[i[0]] = i[2:]

    ############## part2  calculate the similarity of A and B lineage
    sample_num = pd.DataFrame(index = linA_list, columns = linA_list)

    for A in linA_list:
        for B in linA_list:
            if np.isnan(sample_num.loc[B, A]) == True:
                mutat_A = set(Lineage_v[A])
                mutat_B = set(Lineage_v[B])
                at_least_similarity = len(mutat_A & mutat_B) / max(len(mutat_A), len(mutat_B))
                sample_num.loc[B, A] = sample_num.loc[A, B] = at_least_similarity
                
    ############## part3  extract the couple(AB)_"line" simi >=75%
    alpha_list = list(set(["B.1.1.7", "Q.1", "Q.2", "Q.3", "Q.4", "Q.5", "Q.6", "Q.7", "Q.8"]) & set(linA_list))
    simi_list = []
    for a in alpha_list:
        for b in alpha_list:
            if a != b:
                simi_list.append(sample_num.loc[a,b])

    print("\n", "The threathold of cluster is:", str(min(simi_list)), "\n") # 0.7837837837837838

    couple = []
    all_related = []
    for i in sample_num:
        for j in sample_num:
            if j != i and sample_num.loc[i,j] >= min(simi_list):
                couple.append([i,j])
                if i not in all_related: # all_related: all the points who with line(s)
                    all_related.append(i)
                if j not in all_related:
                    all_related.append(j)
            else:
                continue

    ############## part4  for each (has line) point, extract the related points with it.
    group = {}
    for l in all_related:
        each_group = [l]
        for i in couple:
            if i[1] == l:
                each_group.append(i[0])
            elif i[0] == l:
                each_group.append(i[1])
        group[l] = set(each_group)

    threathold_cluster = round(min(simi_list), 4)
    with open(dirpath + "1_filtered_data/simi_group" + str(threathold_cluster) + ".csv","a+") as f:
        for i in group:
            f.write(i + ":")
            for v in group[i]:
                f.write(v + ",")
            f.write("\n")
    ############## part5  for each (has line) point, compare the realted points, whether is the same
    already = []
    concluded_g = {}
    for g in group:
        if g not in already:
            concluded = []
            for m in group:
                if m != g and group[m] == group[g]:
                    concluded.append(m)
                    already.append(m)
                else:
                    continue
            concluded_g[g] = concluded

    ############## part6  for each (has line) point, cheack if they have "close partner(s)", if has, to the same cluster
    count = 0
    cluster = {}
    for i in concluded_g:
        if len(concluded_g[i]) >= 1:
            count += 1
            groupname = "cluster" + str(int(count))
            raw =  list(set(concluded_g[i]))
            raw.append(i)
            raw.sort()
            cluster[groupname] = raw
        else:
            continue

    print("\n", "The total number of cluster is:", count, "\n") # 148

    ############## part7  write the cluster info
    with open(dirpath + "1_filtered_data/cluster" + str(threathold_cluster) + ".csv", "a+") as f:
        for i in cluster:
            f.write(i + ",")
            for v in cluster[i]:
                f.write(v + ",")
            f.write("\n")

    ############## part8  record each points is belong to which cluster
    lin_CLA = {}
    for i in cluster:
        for A in cluster[i]:
            lin_CLA[A] = i

    ############## part9  select the represent points(who with the most patients)
    repres_lin = {}
    for clu in cluster:
        lingeage_n_enrolled = {}
        for lin in linA_list:
            if Lineage_v[lin] >= 1:
                if lin in cluster[clu]:
                    lingeage_n_enrolled[lin] = int(lineage_n[lin])
        repres_lin[clu] = max(lingeage_n_enrolled, key = lingeage_n_enrolled.get)
        
    ############## part10  write down the new FV_75_norm_cluster.txt file
    count_raw = 0
    count_cluster = 0
    cluster_already = []
    for lin in linA_list:
        if lin not in lin_CLA: # individual lineages
            count_raw += 1
            with open(dirpath + "1_filtered_data/"+fv_feature_filename+".txt","a+") as h:
                h.write(lin + "," + lingeage_n_enrolled[lin] + "," + ",".join(Lineage_v[lin]) + "\n")
        else:
            cluster_name = lin_CLA[lin]
            count_cluster += 1
            if cluster_name not in cluster_already:
                cluster_already.append(cluster_name)
                culster_rep_lin = repres_lin[cluster_name]
                culster_num = lineage_n[culster_rep_lin]
                culster_mutations = sort_humanly(set(Lineage_v[culster_rep_lin]))

                with open(dirpath + "1_filtered_data/"+fv_feature_filename, "a+") as h:
                    h.write(culster_rep_lin + "*," + culster_num + "," + ",".join(culster_mutations) + "\n")
    return threathold_cluster


def min_continous_mutations(fv_path, pango_lineage, variants_all):
    from collections import Counter
    from matplotlib.pyplot import plt, MultipleLocator
    # read feature SNP for A pool and B pool
    linA_list = []
    Lineage_v = {}
    with open(fv_path, 'r') as f:
        for i in csv.reader(f): # 2302
            if (i[0].startswith("X")) or (i[0] == "Unassigned") or (i[0] == "None") or (int(i[1]) <= 100):
                continue
            else:
                linA_list.append(i[0])
                Lineage_v[i[0]] = i[2:] # 1491

    # calculate the number of sequential denovo mutations in each sample
    raw_count = {'0': 0, '1': 0}
    for lin in linA_list:
        lin_epi = [k for k, v in pango_lineage.items() if v == lin]
        for i in lin_epi:
            epi_record = []
            num_denovo = 0
            try:
                for v in variants_all[i]:
                    if v in set(Lineage_v[lin]):
                        epi_record.append(num_denovo)
                        num_denovo = 0
                    else:
                        num_denovo += 1
                epi_record.append(num_denovo)
                denovo_count = dict(Counter(epi_record))
                raw_count = dict(Counter(raw_count) + Counter(denovo_count))
            except KeyError:
                continue

    lin_denovo_count = dict([(k, raw_count[k]) for k in sorted(raw_count.keys())])
    # the exsiting chance for different sequential denovo mutations
    num = 0
    for d in lin_denovo_count:
        if d != 0:
            num = num + lin_denovo_count[d]
        else:
            continue
    # the exsiting proporbility for different sequential denovo mutations
    lin_denovo_prop = {}
    for d in lin_denovo_count:
        if d != 0:
            lin_denovo_prop[d] = lin_denovo_count[d] / num
        else:
            continue

    Proportion = list(lin_denovo_prop.keys())[0:13]
    CF = list(lin_denovo_prop.values())[0:13]
    print("------------------- Proportion for different threshold ------------------- " + "\n"+CF+"\n"+Proportion)

    for num in CF:
        if CF[num] < 0.05:
            len_UXY = num

    plt.cla()
    plt.plot(Proportion, CF, color='#8696a7', marker='o')
    plt.xlabel('Number of sequential feature mutations', fontsize=10)
    plt.ylabel('Existing proporbility', fontsize=10)
    ax = plt.gca()
    x_major_locator = MultipleLocator(1)
    ax.xaxis.set_major_locator(x_major_locator)
    y_major_locator = MultipleLocator(0.05)
    ax.yaxis.set_major_locator(y_major_locator)
    plt.grid(True)
    plt.savefig(dirpath + "1_filtered_data/" + "denovo_chance.pdf", bbox_inches='tight')

    return len_UXY
