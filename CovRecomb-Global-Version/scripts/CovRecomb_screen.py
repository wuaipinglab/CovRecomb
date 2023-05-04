import csv
import time
import pandas as pd
import os
import re
from multiprocessing import Process
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests

def recombination_detection(Strain_list_snp, collect_date, variants_all, pango_lineage, Lineage_v, lin_list, lin_earliest_date, must_inA,must_inB,feature_mutations, output_file, mutaions_num):
    dirpath = "/home/wap/jiaying_recom/04_recomb_20230208/"
    os.chdir(dirpath+"scripts/")
    from function_set import dirpath, len_UXY, max_bk_num, bk_count, obtain_pattern,sort_humanly, bk_count,calcul_bk
  
    num = 0
    for epi in Strain_list_snp:
        num += 1
        if num % 1000 == 0 :
            print("Processing {}:    {}/{}    proportion: {}%".format(os.getpid(), num, len(Strain_list_snp), round(num/len(Strain_list_snp),8)*100))
        try:
            epi_time = collect_date[epi]
            epiV = variants_all[epi]
            epi_feat = len(set(epiV) & set(feature_mutations)) # n

            if epi_feat < 2 * len_UXY:
                continue
            else:
                # P-value for Non-recombination
                epi_record = {}
                aftertime_lin = []
                for lin_A in lin_list:
                    timeA = lin_earliest_date[lin_A]
                    if epi_time >= timeA:
                        aftertime_lin.append(lin_A)
                        all_AA = len(Lineage_v[lin_A]) # K
                        all_AA_epi = len(set(Lineage_v[lin_A]) & set(epiV)) # k
                        pVal = hypergeom.sf(all_AA_epi - 1, mutaions_num, all_AA, epi_feat)
                        epi_record[str(lin_A) + "_" + str(lin_A)] = pVal

                # the least p-value for the Non-recombinant
                min_AA = min(epi_record, key = epi_record.get)
                # P-value for Recombinant (A+B/A+B+A)
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
                                            all_AB = len(set(Lineage_v[A]) | set(Lineage_v[B]))  # K
                                            all_AB_epi = len(set(set(Lineage_v[A]) | set(Lineage_v[B])) & set(epiV)) # k
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
                            file_epi.write(epi + "," + collect_date[epi] + "," + pango_lineage[epi]+","+lin_A_draw + "," + lin_B_draw + "," + lin_record + "," +\
                                str(epi_record[min_adjp_pair])+","+str(lin_adjP[min_adjp_pair])+","+"/".join(UA_mutate_unique) + "," + "/".join(UB_mutate_unique) + "," + "/".join(shared_mut) + "," + "/".join(denovo_mut) + "\n")

        except:
            continue

def main():
    dirpath = "/home/wap/jiaying_recom/04_recomb_20230208/"
    os.chdir(dirpath+"scripts/")
    from function_set import dirpath, length, len_UXY, max_bk_num, cor_num, creat_dir
    maintime1 = time.time()
    print("\n", "------------------------- Parameters ------------------------", "\n")

    dirpath_recomb = dirpath + "2_recomb_identified/"
    lineage_file = dirpath + '1_filtered_data/fv_clustered_enrolled.txt'

    print("\n", "The least number of nucleotide in genome: {}".format(length),"\n", \
        "Breakpoint(s) limitation: {}".format(max_bk_num), "\n", \
            "The least number of sequential feature mutations: {}".format(len_UXY), "\n", \
                "The core for running the process: {}".format(cor_num), "\n", \
                    "Start time: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    print("\n", "----------------------- Load datasets -----------------------", "\n")
    output_file = dirpath_recomb + "0_putative_recombinants.csv"
    must_inA, must_inB = "X"*len_UXY, "Y"*len_UXY
    creat_dir(dirpath_recomb)

    ## read the feature mutations and the earliest date for all enrolled (clustered) lineages
    Lineage_v = {}
    lin_list = []
    feature_mutations = []
    num = 0
    with open(lineage_file, 'r') as f:
        for i in csv.reader(f):
            num+=1 
            lin_list.append(i[0])
            Lineage_v[i[0]] = i[2:]
            for v in i[2:]:
                feature_mutations.append(v)

    mutaions_num = len(set(feature_mutations)) # N:
    df = pd.read_csv(dirpath + "1_filtered_data/lineage_earliest_date.csv")
    lin_earliest_date = dict(zip(df['lineage'],df['lineage_earliest_time']))

    ## read the meta file
    meta_path = os.path.join(dirpath, '0_raw_data', 'metadata.tsv')
    with open(meta_path, "r") as f:
        next(f)
        rows = f.readlines()

    num_meta = 0
    pango_lineage = {}
    collect_date = {}
    for i in rows:
        num_meta += 1
        info = i.split("\t")
        seq_length = int(info[8])
        host = info[9]
        if host == "Human" and seq_length >= length and re.search('20\d\d-\d\d-\d\d', info[5]):
            epi = info[0]
            pango_lineage[epi] = info[13]
            collect_date[epi] = info[5]

    del rows
    ## read each sample's mutations
    variants_all = {}
    with open(dirpath + "/1_filtered_data/snp_norm.txt", "r") as f:
        for i in f.readlines():
            i_2 = i.split(':')
            variants_all[i_2[0]] = i_2[1].strip().split(',') 

    Strain_list_snp = list(variants_all.keys())

    print("Input genomes for CovRecomb: ",len(Strain_list_snp),"\n")
    print("Input genomes with metadata: ",len(pango_lineage),"\n")
    print("Candidate parental lineages: ",len(lin_list),"\n") 
    print("Threshold for the continous feature mutations: ",str(len_UXY),"\n") 

    print("\n", "----------------------- CovRecomb analysis -----------------------", "\n")

    col_names = ['sample_id', "collect_date", "pango_lineage", 'lineage_X', 'lineage_Y', \
        'mutation_pattern', "raw_p_value","adjusted_p_value","X_mutations", "Y_mutations", "shared_mutations", "denovo_mutations"]
    with open(output_file, "a+") as file_epi:
        for c in col_names[0:-1]:
            file_epi.write(c+",")
        file_epi.write(col_names[-1]+"\n")

    for i in range(1, cor_num+1):
        globals()["p"+str(i)] = Process(target=recombination_detection, args=(Strain_list_snp[(i-1)*int(len(Strain_list_snp)/cor_num):i*int(len(Strain_list_snp)/cor_num)], \
            collect_date, variants_all, pango_lineage, Lineage_v, lin_list, lin_earliest_date, must_inA, must_inB,\
                feature_mutations, output_file, mutaions_num))
        
    for i in range(1,cor_num+1):
        eval("p"+str(i)).start()
        
    for i in range(1,cor_num+1):
        eval("p"+str(i)).join()
        print("\n","p"+str(i)+" elapsed time： %.8s s" % (time.time() - maintime1),"\n")

    print("\n", "----------------------- Analysis completed -----------------------","\n")


if __name__ == '__main__':
    main()

