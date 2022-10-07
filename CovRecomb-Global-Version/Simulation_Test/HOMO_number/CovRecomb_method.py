import re
from scipy.stats import hypergeom
import copy

def sample_lin_time(lin_A,lin_time):
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


###sort the feature SNP by genome location    
def sort_humanly(v_list): 

    def tryint(s):                       
        try:
            return int(s)
        except ValueError:
            return s

    def str2int(v_str):
        return [tryint(sub_str) for sub_str in re.split('([0-9]+)', v_str)]

    return sorted(v_list, key=str2int,reverse=False)


def calcul_bk_test(lin_A_draw,lin_B_draw,Lineage_v,epiV):
    feature_SNPA = Lineage_v[lin_A_draw]
    feature_SNPB = Lineage_v[lin_B_draw]

    A_B_shared = set(feature_SNPA) & set(feature_SNPB)
    UA_mutate = (set(feature_SNPA) & set(epiV)) - set(A_B_shared)
    UB_mutate = (set(feature_SNPB) & set(epiV)) - set(A_B_shared)

    UA_mutate_unique = []
    UB_mutate_unique = []
    
    lin_record = ""
    for j in epiV:
        if j in UA_mutate:
            UA_mutate_unique.append(j)
            lin_record = lin_record+"A"
        elif j in UB_mutate:
            UB_mutate_unique.append(j)
            lin_record = lin_record+"B"

    return lin_record,UA_mutate_unique,UB_mutate_unique


def unique_lin(aftertime_lin,Lineage_v,U_mutate_unique):
    num = 0
    same_ancient = []
    for can_A in aftertime_lin:
        if set(Lineage_v[can_A]) >= set(U_mutate_unique):
            same_ancient.append(can_A)
            num += 1

    if num ==1:
        return num
    elif num >1:
        same_A_num = mutation_lin_unique(same_ancient)
        if same_A_num == len(same_ancient):
            num =1
        else:
            num =2
    return num


def mutation_lin_unique(same_ancient):
    clean_same_ancient = []
    for mb in same_ancient:
        if "cluster" in mb:
            clean_same_ancient.append(mb.split("_")[1])
        else:
            clean_same_ancient.append(mb)

    same_num = 0
    for saan in clean_same_ancient:
        if len(saan.split("."))<3: 
            break
        else:
            same_01 = clean_same_ancient[0].split(".")[0:3]
            if saan.split(".")[0:3] != same_01:
                break
            elif saan.split(".")[0:3] == same_01:
                same_num +=1
    return same_num


def get_ch1_ch2_name(CH):
    CH1 = CH.split("bk")[0].split("_")[-1]
    if "_" in CH.split("bk")[1]:
        CH2 = CH.split("bk")[1].split("/")[0].split("_")[-1]
    else:
        CH2 = "lin"+CH.split("bk")[1].split("/")[0].split("lin")[-1]
    return CH1,CH2


def get_lin_name(CH):
    if "gen" in CH:
        CH1 = CH2 = "lin" + CH.split("_lin")[1]
    else:
        CH1 = CH2 = CH
    return CH1,CH2


#################
def recombination_detection_test(Strain_list_snp,len_UAB,linA_list,variants_all,feature_mutations,Lineage_v,mutaions_num):
    must_inA = "A" * len_UAB
    must_inB = "B" * len_UAB
    linA_list_deep = copy.deepcopy(linA_list)
        
    correct_AA_num = 0
    correct_AB_num = 0
    false_AA_num = 0
    false_AB_num = 0

    zhenyin = 0
    jiayang = 0
    jiayin = 0
    zhenyang = 0
    for epi in Strain_list_snp:
        if "bk" in epi:
            CH1,CH2 = get_ch1_ch2_name(epi)
        else:
            CH1,CH2 = get_lin_name(epi)

        try:
            epi_record = {}
            epiV = variants_all[epi]
            epi_feat = len(set(epiV) & set(feature_mutations))

            ### P-value for Non-recombination
            for lin_A in linA_list_deep:
                all_AA = len(Lineage_v[lin_A])
                all_AA_epi = len(set(Lineage_v[lin_A]) & set(epiV))

                pVal = hypergeom.logsf(all_AA_epi-1,mutaions_num,all_AA,epi_feat)
                epi_record[str(lin_A)+"_"+str(lin_A)] = (pVal*10000000000)

            # the least p-value for the Non-recombinant
            most_one = min_pairs(epi_record)
            epi_record  = {}
            for mo in most_one:
                epi_record[mo[0]] = mo[1]

            ### P-value for Recombinant (A+B/A+B+A)
            A_already = []
            for A in linA_list_deep:
                A_already.append(A)
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
   
                                epi_record[str(A)+"_"+str(B)] = (pVal*10000000000)
                
            epi_record_sort = sorted(epi_record.items(), key=lambda item:item[1], reverse=False)
        except:
            continue
        
        most_two = min_pairs(epi_record)
        
        if (CH1 == CH2) and (most_one == most_two): # Not recombinant or intra-lineage recombinant, recognized as not recombinant (CORRECT)
            zhenyin += 1
            if {most_two[0][0].split("_")[0],most_two[0][0].split("_")[1]} - {CH1,CH2} == set(): # Not recombinant, recognized as not recombinant (CORRECT AA)
                correct_AA_num += 1
                continue
            else:
                false_AA_num += 1 # Not recombinant, recognized as not recombinant (FALSE AA)
                continue
            
        elif CH1 == CH2 and (most_one != most_two): # Not recombinant, but recognized as recombinant (FALSE)
            jiayang += 1
            false_AA_num +=1

        elif CH1 != CH2 and (most_one == most_two):  # Recombinant, but recognized as not recombinant (FALSE)
            jiayin += 1
            false_AB_num += 1 

        elif CH1 != CH2 and (most_one != most_two):  # Recombinant, recognized as recombinant (CORRECT)
            zhenyang += 1
            most_two = [epi_record_sort[0]] # AB
            result_recom = {str(most_two[0][0].split("_")[0]), str(most_two[0][0].split("_")[1])}
            epiV = sort_humanly(epiV)

            select_recom = most_two[0][0].split("_")
            lin_A_draw = select_recom[0]
            lin_B_draw = select_recom[1]

            lin_record, UA_mutate_unique, UB_mutate_unique = calcul_bk_test(lin_A_draw,lin_B_draw,Lineage_v,epiV)
            numA = unique_lin(linA_list_deep, Lineage_v, UA_mutate_unique)
            numB = unique_lin(linA_list_deep, Lineage_v, UB_mutate_unique)
            if numA == numB == 1:
                if result_recom - {CH1, CH2} == set():
                    correct_AB_num += 1 # Recombinant, recognized as recombinant (CORRECT AB)
            else:
                false_AB_num += 1 # Recombinant, recognized as recombinant (FALSE AB)
        
    sample_num = len(Strain_list_snp)
    return sample_num, correct_AB_num, false_AB_num, correct_AA_num, false_AA_num, zhenyin, jiayang, jiayin, zhenyang
