'''
Function: Use the CovRecomb method to identify inter-lineage recombinants.

Input: Four files: snp_norm.txt, fv_norm_cluster.txt, Lineage_Earliest_date.csv, metadata-2.tsv, 

Output: One file: 0_putative_recombinants.csv, record the idenfied putative recombinants with their parental lineages and feature mutation patterns.
'''

import csv
import time
import pandas as pd
import os
import argparse
from multiprocessing import Process

def main():
    parser = argparse.ArgumentParser(description='Read the qulified sequences.', usage='''python3 CovRecomb_pipeline.py''')
    
    parser.add_argument("-dir", "--dirpath", default="/home/soniali/Desktop/03_CovRecomb_review0819/github/CovRecomb/CovRecomb-Global-Version/demo_for_hypergeometric/data/2022_02_12/", help="The totoal folder address of this work.\nDefault: /home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12")
    parser.add_argument("-b", "--breakpoint_number", default=2, help="The maximum acceptable number of breakpoints among putative recombinants.\nDefault: 2")
    parser.add_argument("-t", "--threshold", default=4, help ="The least number of sequential feature mutations.\nDefault: 4")
    parser.add_argument("-c", "--core_number", default=16, help="The number of cores used while computation running.\nDefault: 4")
    parser.add_argument("-o", "--output_file", default="0_putative_recombinants.csv", help="The filename of the output file with identified putative recombinants and their parental sequences.\nDefault: putative_recombinants.csv\n")
    args = parser.parse_args()
    
    dirpath = args.dirpath
    # dirpath = "/home/soniali/Desktop/03_CovRecomb_review0819/github/CovRecomb-Global-Version/demo/data/2022_02_12/"
    os.chdir(dirpath.split("data/2022_02_12/")[0]+"scripts/")
    from functions_set import creat_dir
    from CovRecomb_detecion import recombination_detection

    maintime = time.time()
    dirpath_recomb = dirpath + "2_recomb_identified/"
    lineage_file = dirpath + '1_filtered_data/fv_norm_cluster.txt'
    output_file = dirpath_recomb + args.output_file
    # output_file = dirpath_recomb + "0_putative_recombinants.csv"
    max_bk_num = int(args.breakpoint_number)
    len_UAB = int(args.threshold)
    cor_num = int(args.core_number)
    
    creat_dir(dirpath_recomb)
    print("\n", "----------------------- Load datasets -----------------------", "\n")
    
    Lineage_v = {}
    linA_list = []
    feature_mutations = []
    with open(lineage_file, 'r') as f:
        for i in csv.reader(f):
            linA_list.append(i[0])
            Lineage_v[i[0]] = i[2:]
            for v in i[2:]:
                if v not in feature_mutations:
                    feature_mutations.append(v)

    mutaions_num = len(feature_mutations)
  
    lin_time = {}
    with open(dirpath + "0_raw_data/Lineage_Earliest_date.csv", "r") as f_lintime:
        f_csv = csv.reader(f_lintime)
        next(f_csv)
        for row in f_csv:
            lin_time[row[0]] = row[1]
            
    # read each sample's mutations
    Strain_list_snp = []
    variants_all = {}
    with open(dirpath + "1_filtered_data/snp_norm_demo.txt", "r") as f:
        for i in f.readlines():
            i_2 = i.split(':')
            Strain_list_snp.append(i_2[0])
            variants_all[i_2[0]] = i_2[1].strip().split(',')
            
    # read each sample's epi_info
    Strain_list_meta = []
    collected_time = {}
    region = {}
    country = {}
    division = {} 
    pango_lineage = {}
    with open(dirpath + "0_raw_data/metadata-demo.tsv", "r") as f:
        next(f)
        for i in csv.reader(f, delimiter = '\t'):
            Strain_list_meta.append(i[2])
            collected_time[i[2]] = i[4]
            region[i[2]] = i[9]
            country[i[2]] = i[10]
            division[i[2]] = i[11]
            pango_lineage[i[2]] = i[18]

    # Strain_list_snp = Strain_list_snp[int(len(Strain_list_snp)/1000*560)+1:int(len(Strain_list_snp)/1000*581)+1]
    
    must_inA = "X"*len_UAB
    must_inB = "Y"*len_UAB
    
    print("Input genomes for CovRecomb: ",len(Strain_list_snp),"\n")
    print("Input genomes with metadata: ",len(Strain_list_meta),"\n")
    print("Candidate parental lineages: ",len(linA_list),"\n")
    print("Threshold for the continous feature mutations: ",str(len_UAB),"\n")
    
    print("\n", "----------------------- CovRecomb analysis -----------------------", "\n")
    
    col_names = ['sample_id', 'lineage_X', 'lineage_Y', 'mutation_pattern', "X_mutations", "Y_mutations", "shared_mutations", "denovo_mutations"]
    with open(output_file, "a+") as file_epi:
        for c in col_names[0:-1]:
            file_epi.write(c+",")
        file_epi.write(col_names[-1]+"\n")

    for i in range(1,cor_num+1):
        globals()["p"+str(i)] = Process(target=recombination_detection, args=(len_UAB, max_bk_num, must_inA, must_inB, linA_list,\
            Strain_list_snp[(i-1)*int(len(Strain_list_snp)/cor_num):i*int(len(Strain_list_snp)/cor_num)],\
            collected_time, variants_all, feature_mutations, lin_time, Lineage_v, mutaions_num, output_file))
        
    for i in range(1,cor_num+1):
        eval("p"+str(i)).start()
        
    for i in range(1,cor_num+1):
        eval("p"+str(i)).join()
        print("\n","p"+str(i)+ " elapsed time： %.8s s" % (time.time() - maintime),"\n")

    print("\n","----------------------- Analysis completed -----------------------","\n")
    
    if os.path.exists(output_file):
        df = pd.read_csv(output_file)
        print("The number of putative inter-lineage recombinants: ", str(df.shape[0]) ,"\n")
    else:
        print("The number of putative inter-lineage recombinants: 0", "\n")
        
    print("Total elapsed time：%.8s s" % (time.time() - maintime),"\n")
        

if __name__ == '__main__':
    main()
