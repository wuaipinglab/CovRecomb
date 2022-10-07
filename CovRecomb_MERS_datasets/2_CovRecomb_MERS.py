
import csv
import re
import pandas as pd
from scipy.stats import hypergeom
import copy

####
def mult_times(str, n):
    must_inA = ""
    for i in range(n):
        must_inA = must_inA+str
    return must_inA

def min_pairs(dic):
    if len(dic) == 0:
        return []
    min_val = min(map(lambda v: v[1], dic.items()))
    return [item for item in dic.items() if item[1] == min_val]

###sort the feature SNP by genome location
def tryint(s):                       
    try:
        return int(s)
    except ValueError:
        return s


def str2int(v_str):     
    return [tryint(sub_str) for sub_str in re.split('([0-9]+)', v_str)]
    
    
def sort_humanly(v_list):   
    return sorted(v_list, key=str2int,reverse=False)


infile = "/home/soniali/Desktop/03_CovRecomb/CovRecomb/MERS/"
DIRPATH = infile
Strain_list_snp = []
variants_all = {}
with open(DIRPATH+"snp.txt","r") as f:
    for i in f.readlines(): 
        i_2 = i.split(':')  
        Strain_list_snp.append(i_2[0])
        variants_all[i_2[0]] = i_2[1].strip().split(',')[:-1]

Strain_list_meta = []
pango_lineage = {}
with open(infile+"group.csv","r") as f:
    next(f)
    for i in csv.reader(f,delimiter=','):
        Strain_list_meta.append(i[2])
        pango_lineage[i[2]] = i[1]

for l in [k for k in pango_lineage.keys() if ".1" in k or ".2" in k]:
    l_short = l.split(".")[0]
    pango_lineage[l_short] = pango_lineage[l]
    
lineage_file = 'fv_75.txt'
Lineage_v = {}
linA_list = []
feature_mutations = []
with open(DIRPATH+lineage_file,'r') as f:
    count = 0
    for i in csv.reader(f):
        if i[2] != "":
            if len(i[2:])>= 1:
                linA_list.append(i[0])
                Lineage_v[i[0]] = i[2:-1]
                for v in i[2:]:
                    if v not in feature_mutations:
                        feature_mutations.append(v)

mutaions_num = len(feature_mutations)
len_UAB = int(4)
must_inA = mult_times("A",len_UAB)
must_inB = mult_times("B",len_UAB)
linA_list_deep = copy.deepcopy(linA_list)

epi_need = []
for epi in Strain_list_snp:
    # epi = "KT368889.1"
    epi_record = {}
    epiV = variants_all[epi]
    epi_feat =len(set(epiV) & set(feature_mutations))
    ### P-value for Non-recombination
    for lin_A in linA_list_deep:
        all_AA = len(Lineage_v[lin_A])
        all_AA_epi = len(set(Lineage_v[lin_A]) & set(epiV))

        pVal = hypergeom.logsf(all_AA_epi-1,mutaions_num,all_AA,epi_feat)
        epi_record[str(lin_A)+"_"+str(lin_A)] = pVal

    # the least p-value for the Non-recombinant
    most_one = min_pairs(epi_record)
    pmin = most_one[0][1]
    epi_record  = {}
    for mo in most_one:
        epi_record[mo[0]] = mo[1]

    ### P-value for Recombinant (A+B/A+B+A)
    A_already = []
    for A in linA_list_deep:
        A_epi = set(Lineage_v[A]) & set(epiV)
        if len(A_epi) < len_UAB:
            continue
        else:
            aftertime_linB = set(linA_list_deep) - set(A_already)
            for B in aftertime_linB:
                B_epi = set(Lineage_v[B]) & set(epiV)
                if len(B_epi) < len_UAB:
                    continue
                elif len(B_epi - A_epi) < len_UAB or len(A_epi-B_epi) < len_UAB:
                    continue
                else:
                    all_AB = (len(Lineage_v[A]) + len(Lineage_v[B]))/2
                    all_AB_epi = len((set(Lineage_v[A]) | set(Lineage_v[B])) & set(epiV))
                    pVal = hypergeom.logsf(all_AB_epi-1,mutaions_num,all_AB,epi_feat)
                    union_AB_set = set(Lineage_v[A]) ^ set(Lineage_v[B])
                    unique_A = A_epi-B_epi
                    unique_B = B_epi-A_epi
                    AB_epi =sort_humanly(list(union_AB_set &  set(epiV)))
                    
                    recom_model = ""
                    for v in AB_epi:
                        if v in unique_A:
                            recom_model = recom_model+"A"
                        elif v in unique_B:
                            recom_model = recom_model+"B"

                    if (must_inA not in recom_model) or (must_inB not in recom_model):
                        continue
                    else:
                        start = recom_model[0]
                        change = 0
                        for R in recom_model:
                            if R != start:
                                change +=1
                                start = R

                        epi_record[str(A)+"_"+str(B)] = pVal
    
    epi_record_sort = sorted(epi_record.items(), key=lambda item:item[1], reverse=False)
    most_two = min_pairs(epi_record)
    if most_one[0] == most_two[0]: # Mostly, this sample is non-recombinant
        continue
    else:
        epi_need.append(epi)
        if ".1" in epi or ".2" in epi:
            epi = epi.split(".")[0]
            print(epi,":",pango_lineage[epi],"_",epi_record_sort)
