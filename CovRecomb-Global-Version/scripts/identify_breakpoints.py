'''
Function: Detect the breakpoint distributions of independent recombination events by the help of the 3SEQ method.

Input: Four files: Independent_recombination_event_94.csv

Output: Three files: (1) seq_results_.txt; the breakpoint distribution identified by 3SEQ
				     (2) trimed_94recom.fasta; the aligned and trimed fasta file for independent recombination events
                     (3) parental_.fasta, child_.fasta; two artificial parental lineages' sequences composed by feature mutations and the real sequence for one recombinant.
'''

import pandas as pd
from Bio import SeqIO
import pandas as pd
import os
import subprocess
import argparse

def get_genome_sites(sequence):
    sites = []
    n = 0
    for i in sequence:
        if i != '-':
            sites.append(n)
        n += 1
    return sites


def trim_seq(sequence, sites):
    seq_trimed_list = []
    for i in sites:
        seq_trimed_list.append(str.upper(sequence[i]))
    seq_trimed = ''.join(seq_trimed_list)
    return seq_trimed


def creat_dir(new_dir):
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
            
            
def main():
    parser = argparse.ArgumentParser(description='Read the qulified sequences.', usage='''python3 identify_patterns.py''')

    parser.add_argument("-dir", "--dirpath", default="/home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12/", help="The totoal folder address of this work.\nDefault: /home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12")
    parser.add_argument("-ev", "--events_file", default="Independent_recombination_event_94.csv", help="The filename of the putative recombinants already verified by the epidemiology background.\nDefault: 1_putative_recomb_final.csv\n")
    args = parser.parse_args()

    dirpath = args.dirpath
    lineage_file = dirpath + '1_filtered_data/fv_norm_cluster.txt'
    events_file = args.events_file
    os.chdir(dirpath.split("data/2022_02_12/")[0] + "scripts/")

    ### Read the independent recombination events
    dirpath_pattern = dirpath + "3_recom_pattern/"
    outpath = dirpath_pattern+"BK_3SEQ/"
    creat_dir(outpath)
    df = pd.read_csv(dirpath_pattern+events_file)

    recom_id_AB = {}
    all_AB = []
    for i in df.index:
        epi = df.loc[i, "sample_id"]
        A_lin, B_lin = df.loc[i, "lineage_X"], df.loc[i, "lineage_Y"]
        recom_id_AB[epi] = [A_lin, B_lin]
        all_AB.append(A_lin)
        all_AB.append(B_lin)

    all_AB = list(set(all_AB))

    # step1 : Acquire the fasta file for events
    df_meta = pd.read_csv(outpath+"Recom_94_meta.csv")
    trimed_file_name = "trimed_94recom"
    rename_sample = {}
    for i in df_meta.index:
        rename_sample[df_meta.loc[i,"Virus name"]] = df_meta.loc[i,"Accession ID"]

    with open(outpath+"maffted_94recom.fasta", "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            ID_name = str(record.id)
            ID_name = ID_name.split("|EPI_ISL")[0]
            if ID_name in rename_sample:
                seq = str(record.seq)
                seq_upper = str(seq)
                rename_id = rename_sample[ID_name]
                with open(outpath+"maffted_94recom_epi.fasta", "a+") as h:
                    h.write('>'+rename_id+"\n")
                    h.write(seq_upper+"\n")

    # step2 : get the reference genomes' sites
    with open(outpath+"maffted_94recom_epi.fasta", "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            ID_name = record.id
            if ID_name == "EPI_ISL_402125":
                seq = str(record.seq)
                sites = get_genome_sites(seq)
                break

    # step3 : Trim the gisaid download fasta sequences
    num = 0
    within_fasta_epi = []
    with open(outpath+"maffted_94recom_epi.fasta", "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            ID_name = str(record.id)
            within_fasta_epi.append(ID_name)
            seq = str(record.seq)
            seq_trimed = trim_seq(seq, sites)
            with open(outpath+trimed_file_name+".fasta", "a+") as fw:
                num += 1
                fw.write(">"+ID_name+"\n")
                fw.write(seq_trimed+"\n")

    ref_genome = []
    with open(outpath+trimed_file_name+".fasta") as f:
        for record in SeqIO.parse(f, "fasta"):
            ID_name = str(record.id)
            if ID_name == "EPI_ISL_402125":
                seq = str(record.seq)
                ref_genome = str.upper(seq)

    # step4 : Generate the artifical sequences composed by lineages' feature mutations
    lineage_mut = {}
    with open(lineage_file) as f:
        lines = f.readlines()
        for i in lines:
            lineage = i.strip().split(",")[0]
            if "cluster" in lineage:
                lineage = lineage.split("_")[1]+"*"
            if lineage in all_AB:
                lineage_mut[lineage] = i.strip().split(",")[2:]

    lineage_mut_site = {}
    lineage_mut_site = {}
    for lin in lineage_mut:
        mut_site = []
        for site in lineage_mut[lin]:
            mut_site.append(int(site.split("_")[0]))
        len(mut_site)
        lineage_mut_site[lin] = mut_site

    lin_seq_half = {}
    for lin in lineage_mut:
        lin_seq = ""
        jump = 0
        need_jump = 0
        for pos in range(len(ref_genome)):
            if jump == need_jump:
                if pos+1 in lineage_mut_site[lin]:
                    for i in lineage_mut[lin]:
                        if str(pos+1) == i.split("_")[0] and "---" not in i:
                            letter = i.split('_')[1]
                            lin_seq += letter
                        elif str(pos+1) == i.split("_")[0] and "---" in i:
                            need_jump = len(i.split('_')[1]) - 1
                            letter = i.split('_')[1]
                            lin_seq += letter
                else:
                    lin_seq += ref_genome[pos]
            else:
                need_jump = need_jump - 1
                continue
        lin_seq_half[lin] = lin_seq
    
    # step5 : Write the triplet sequences (two parental lineages + one child sequence)
    included_name = []
    with open(outpath+trimed_file_name+".fasta") as f:
        for record in SeqIO.parse(f, "fasta"):
            ID_name = str(record.id)
            if ID_name in recom_id_AB:
                seq_epi = str(record.seq)
                epi_path = outpath+ID_name+"/"
                creat_dir(epi_path)
                included_name.append(ID_name)
                with open(epi_path+ID_name+"_child.fasta","a+") as f:
                    f.write(">"+ID_name+"\n")
                    f.write(seq_epi+"\n")
                    
                with open(epi_path+ID_name+"_parental.fasta","a+") as f:
                    for lin in recom_id_AB[ID_name]:
                        f.write(">"+lin+"\n")
                        f.write(lin_seq_half[lin]+"\n")

    # step6 : With the help of 3SEQ to identify breakpoint distribution
    for epi in recom_id_AB:
        path = outpath+epi+"/"
        if os.path.exists(path) and len(os.listdir(path)) == 2:
            os.chdir(path)
            for fi in os.listdir(path):
                if "_parental.fasta" in fi:
                    parene_file_path = path+"/"+fi
                    record_para = path+"/"+fi.split("_parental_seq")[0]
                elif "_child.fasta" in fi:
                    child_file_path = path+"/"+fi
                    record_para = path+"/"+fi.split("_recom_seq")[0]

            subprocess.call (["cd /home/soniali/Downloads/3seq_build/ && \
                echo y | ./3seq -f %s %s -id %s -t1" % (parene_file_path, child_file_path,record_para)],shell=True)

    # step7 : Summarize the 3SEQ resutls:
    seq_results = []
    seq_without = []
    for epi in recom_id_AB:
        path = outpath+epi+"/" 
        for file in os.listdir(path):
            if "fasta.3s.rec" in file:
                with open(path+file) as f:
                    lines = f.readlines()
                    if len(lines) == 2:
                        seq_results.append(lines[1])
                    else:
                        seq_without.append(epi)
                        
    print("\n", "----------- Sequence(s) failed to be identified as recombinants by 3SEQ ----------------", "\n")
    print(seq_without, "\n")
    
    print("\n", "----------- Write the breakpoint distribution identified by 3SEQ ----------------", "\n")
    with open(outpath+"seq_results_"+str(len(seq_results))+".txt", "a+") as f:
        for i in seq_results:
            f.write(i)


if __name__ == '__main__':
    main()
    