'''
The collection set of functions
'''
import pandas as pd
import copy
import re
import os
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import hypergeom

# length = 27000
# fm_threshold = 0.75
# len_UXY = 4
# cor_num = 20
# dirpath = "/home/soniali/Desktop/02_recom_230203/data/2023_02_08/"
# sequence_path = os.path.join(dirpath, '0_raw_data', 'nextclade.tsv')
# meta_path = os.path.join(dirpath, '0_raw_data', 'metadata.tsv')
# geno_min,geno_max = 266,29674
# mini_simi_default = 0.7838
# snps_path = dirpath+"/1_filtered_data/snp_norm.txt"
# fv_feature_filename = "fv_clustered_enrolled.txt"

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


'''
The core algorithm for CovRecomb method
'''
def recombination_detection(Self_snp, Self_variants, Lineage_v, lin_list, must_inA,must_inB,feature_mutations, output_file, mutaions_num,max_bk_num,len_UXY):
    for epi in Self_snp:
        epi
        try:
            epiV = Self_variants[epi]
            epi_feat = len(set(epiV) & set(feature_mutations))
            # P-value for Non-recombination
            epi_record = {}
            aftertime_lin = []
            for lin_A in lin_list:
                aftertime_lin.append(lin_A)
                all_AA = len(Lineage_v[lin_A]) 
                all_AA_epi = len(set(Lineage_v[lin_A]) & set(epiV))
                pVal = hypergeom.sf(all_AA_epi - 1, mutaions_num, all_AA, epi_feat)
                epi_record[str(lin_A) + "_" + str(lin_A)] = pVal

            # the least p-value for the Non-recombinant
            min_AA = min(epi_record, key = epi_record.get)
            A_already = []
            for A in aftertime_lin:
                A_already.append(A)
                A_epi = set(Lineage_v[A]) & set(epiV)
                if len(A_epi) < len_UXY:
                    continue
                else:
                    afterA_linB = set(aftertime_lin) - set(A_already)
                    for B in afterA_linB:
                        B_epi = set(Lineage_v[B]) & set(epiV)
                        if len(B_epi) < len_UXY:
                            continue
                        else:
                            unique_A = A_epi - B_epi
                            unique_B = B_epi - A_epi
                            if len(unique_A) < len_UXY or len(unique_B) < len_UXY:
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
                                        all_AB = len(set(Lineage_v[A]) | set(Lineage_v[B]))  
                                        all_AB_epi = len(set(set(Lineage_v[A]) | set(Lineage_v[B])) & set(epiV)) 
                                        pVal = hypergeom.sf(all_AB_epi - 1, mutaions_num, all_AB, epi_feat)
                                        epi_record[str(A) + "_" + str(B)] = pVal

            raw_pvals = list(epi_record.values())
            rejected, p_adjusted, _, alpha_corrected = multipletests(raw_pvals, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)

            lin_adjP = {}
            for p in range(len(p_adjusted)):
                lin_pair = list(epi_record.keys())[p]
                lin_adjP[lin_pair] = p_adjusted[p]
            
            min_adjp_pair = min(lin_adjP, key = lin_adjP.get)
            if min_adjp_pair == min_AA or lin_adjP[min_adjp_pair] >= 0.05: 
                continue
            else:
                epiV = sort_humanly(epiV)
                lin_A_draw, lin_B_draw = min_adjp_pair.split("_")[0],min_adjp_pair.split("_")[1]
                lin_record, UA_mutate_unique, UB_mutate_unique, shared_mut, denovo_mut = calcul_bk(lin_A_draw, lin_B_draw, Lineage_v, epiV)
                if set(UA_mutate_unique + UB_mutate_unique) - set(Lineage_v[min_AA.split("_")[0]]) == set():
                    continue
                else:
                    with open(output_file, "a+") as file_epi:
                        file_epi.write(epi + "," +lin_A_draw + "," + lin_B_draw + "," + lin_record + "," +\
                            str(epi_record[min_adjp_pair])+","+str(lin_adjP[min_adjp_pair])+","+"/".join(UA_mutate_unique) + "," + "/".join(UB_mutate_unique) + "," + "/".join(shared_mut) + "," + "/".join(denovo_mut) + "\n")
        except:
            continue