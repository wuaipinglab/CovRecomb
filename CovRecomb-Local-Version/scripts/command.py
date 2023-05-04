'''
The main command for CovRecomb-Local-Version to detect recombinant(s) among all inputted sequences.
'''

import pandas as pd
import csv
import time
import os
import warnings
from Bio import SeqIO
import argparse
from multiprocessing import Pool, Process

os.chdir(os.getcwd() + "/")
from function_set import normalize_position, recombination_detection

def read_seq(sequence_path, REF, times):
    '''
    Read the sequences in the input .fasta file.

    :param sequence_path: str, the address of the .fasta file
    :param date_2: REF, the name of the reference genome in the .fasta file
    :param times: int, the times number of the read_seq process

    :return: boolean, whether the process of reading all the genomes in the .fasta file has been finished
    :return: dict, key: name of each sequence; value: mutations
    '''
    # print('Read times: ', times)

    msa_seqs = []
    seqs = {}
    with open(sequence_path) as f:
        for n, seq_record in enumerate(SeqIO.parse(f, 'fasta')):
            seq_id = seq_record.description
            nstart = 500000 * (times - 1)
            nend = 500000 * times
            if nstart <= n < nend or seq_id == REF:
                msa_seqs.append(seq_id)
                seqs[seq_id] = str(seq_record.seq).upper()
    print("\n", 'Number of input genomes: ', len(seqs), "\n")

    if nend > n:
        finish_all_seq = True
    else:
        finish_all_seq = False

    return finish_all_seq, seqs


def get_genome_sites(sequence):
    '''
    Get the genome position for each site in sequence

    :param sequence: str, the genome sequence composed of 'A, T, C, G, -'

    :return: list, all positions with 'A, T, C, G' bases
    '''
    sites = []
    n = 0
    for i in sequence:
        if i != '-':
            sites.append(n)
        n += 1
    return sites


def trim_seq(sequence, sites):
    '''
    Trim the genome sequence

    :param sequence: str, the genome sequence need to be trimed
    :param sites: int, the genome position need to be preserved

    :return: list, the trimed sequence
    '''
    seq_trimed_list = []
    for i in sites:
        seq_trimed_list.append(sequence[i])
    seq_trimed = ''.join(seq_trimed_list)
    return seq_trimed


def is_qualified(sequence):
    '''
    Check whether the sequence is qualified: 'ATCG >= 27000'

    :param sequence: str, the genome sequence need to be checked

    :return: boolean, whether the sequence is qualified
    '''
    if sequence.count('A') + sequence.count('T') + sequence.count(
            'G') + sequence.count('C') >= 27000:
        return True
    else:
        return False


def get_snps(sequence, ref_genome):
    '''
    Get the mutations for each qualified sequence

    :param sequence: str, the genome sequence needs to be analyzed
    :param ref_genome: str, the genome of the reference sequence 

    :return: list, the mutations of the sequence
    '''
    snps = []
    gaps = []
    for gsite in range(0, len(ref_genome)):
        if sequence[gsite] != ref_genome[gsite] and sequence[gsite] in [
                'A', 'T', 'G', 'C'
        ]:
            snps.append(str(gsite + 1) + '_' + sequence[gsite])
        elif sequence[gsite] == '-':
            gaps.append(gsite + 1)

    continuous_gaps = []
    for gap in gaps:
        if len(continuous_gaps) == 0 or gap == continuous_gaps[-1] + 1:
            continuous_gaps.append(gap)
        else:
            if len(continuous_gaps) % 3 == 0 and continuous_gaps[
                    0] != 1 and continuous_gaps[-1] != 29891:
                snps.append(
                    str(continuous_gaps[0]) + '_' + '-' * len(continuous_gaps))
            continuous_gaps = [gap]

    return snps


def clean_seq(seq_id, sequence, sites, ref_genome):
    '''
    Filter out the qualified sequence and get the mutaions

    :param seq_id: str, the name of each sequence in the input .fasta file
    :param sequence: str, the genome sequence needs to be analyzed
    :param sites: str, the trimed sites of the reference sequence 
    :param ref_genome: str, the genome of the reference sequence 

    :return: dict, key: the name of each sequence; value: the mutations of each sequence
    '''
    seq_trimed = trim_seq(sequence, sites)
    if is_qualified(seq_trimed):
        snps = get_snps(seq_trimed, ref_genome)
        return {seq_id: snps}
    return False


def feature_mut(DIRPATH, lineage_file):
    '''
    Get the feature mutations for each representative lineage

    :param DIRPATH: str, the file address of the file below
    :param lineage_file: str, the filename of file which contains the lineage-defining feature mutations (LDFM) library

    :return: list, The candidate parental lineages
    :return: dict, key: the candidate parental lineage; value: the feature mutations
    :return: int, the number of all the unique mutations in the LDFM library
    :return: list, all the unique mutations in the LDFM library
    '''

    Lineage_v = {}
    lin_list = []
    feature_mutations = []
    with open(DIRPATH + lineage_file, 'r') as f:
        for i in csv.reader(f):
            if i[2] != "":
                if len(i[2:]) >= 1:
                    lin_list.append(i[0])
                    Lineage_v[i[0]] = i[2:]
                    for v in i[2:]:
                        if v not in feature_mutations:
                            feature_mutations.append(v)

    mutaions_num = len(feature_mutations)

    return lin_list, Lineage_v, mutaions_num, feature_mutations


def mafft_to_snp_norm(sequence_path,
                      self_path,
                      len_UXY,
                      cor_num,
                      REF="EPI_ISL_402125"):
    '''
    Get the normalized mutations for sequences in the input .fasta file

    :param sequence_path: str, the file address and the filename of the input .fasta file
    :param self_path: the file address of the input .fasta file
    :param len_UXY: int, the least number of the sequential feature mutations  
    :param REF: str, the name of the reference seq
    '''
    snps_norm_path = self_path + "snp_norm.txt"
    warnings.filterwarnings('ignore')

    times = 0
    finish_all_seq = False
    seq_snps = {}
    while not finish_all_seq:
        times += 1

        finish_all_seq, seqs = read_seq(sequence_path, REF, times)
        sites = get_genome_sites(seqs[REF])
        ref_genome = trim_seq(seqs[REF], sites)

        result_list = []
        p = Pool(cor_num)
        for i in seqs:
            if i != REF:
                result = p.apply_async(clean_seq,
                                       args=(i, seqs[i], sites, ref_genome))
                result_list.append(result)
        p.close()
        p.join()
        for result in result_list:
            if result.get():
                seq_snps.update(result.get())

    print("\n", 'The number of genomes with high quality: ', len(seq_snps),
          "\n")
    for id in seq_snps:
        mut = seq_snps[id]
        nor_snp = normalize_position(mut, 266, 29674)
        if nor_snp.count(",") + 1 < (2 * len_UXY):
            continue
        else:
            with open(snps_norm_path, "a+") as fsn:
                fsn.write(id + ":" +
                          "".join(nor_snp[1:-1].replace("'", "").split()) +
                          "\n")


from scripts import __version__


def main():
    # command line interface
    parser = argparse.ArgumentParser(
        description='The local version of CovRecomb method to detect inter-lineage recombinants.',
        usage='''CovRecomb <your_address_and_filename.fasta>''')
    parser.add_argument(
        'alignment',
        help="The filename of the input fasta file with aligned sequences, including the reference sequence."
    )
    parser.add_argument(
        "-f",
        "--file_address",
        default=os.getcwd() + "/",
        help="The address of the input file.\nDefault: /CovRecomb-Local-Version/")
    parser.add_argument(
        "-b",
        "--breakpoint_number",
        default=2,
        help="The maximum acceptable number of breakpoints among putative recombinants.\nDefault: 2"
    )
    parser.add_argument(
        "-r",
        "--reference",
        default="EPI_ISL_402125",
        help="The reference sequence ID, indicating which sequence\nin the alignment is the reference (by sequence ID).\nDefault: EPI_ISL_402125"
    )
    parser.add_argument(
        "-t",
        "--threshold",
        default=4,
        help="The least number of sequential feature mutations.\nDefault: 4")
    parser.add_argument(
        "-o",
        "--output_file",
        default="putative_recombinants.csv",
        help="The filename of the output file with identified putative recombinants and their parental sequences.\nDefault: putative_recombinants.csv\n"
    )
    parser.add_argument(
        "-c",
        "--core_number",
        default=2,
        help="The number of cores used while computation running.\nDefault: 2")
    args = parser.parse_args()

    maintime = time.time()
    download_path = args.file_address  #"example_project/input_file/Example_aligned_trimed.fasta"
    file_path = args.alignment
    sequence_path = download_path + file_path
    output_file = download_path + args.output_file  

    REF = args.reference
    max_bk_num = int(args.breakpoint_number)
    len_UXY = int(args.threshold)
    cor_num = int(args.core_number)

    print("\n", "----------------------- STEP1: Date filtration and mutations extraction ------------------")
    # get the address of the input file
    self_path = download_path
    for i in range(len(file_path.split("/")) - 1):
        self_path += file_path.split("/")[i] + "/"

    # get the mutations for each sequence
    if os.path.exists(self_path + "snp_norm.txt"):
        os.remove(self_path + "snp_norm.txt")

    mafft_to_snp_norm(sequence_path, self_path, len_UXY, cor_num, REF)
    Self_snp = []
    Self_variants = {}
    with open(self_path + "snp_norm.txt", "r") as f:
        for i in f.readlines():
            i_2 = i.split(':')
            Self_snp.append(i_2[0])
            Self_variants[i_2[0]] = i_2[1].strip().split(',')

    print("\n", "----------------------- STEP2: Get the lineage-defining feature mutations -------------", "\n")
    lineage_file = 'defaults/LDFM_Feb08_2023.txt'
    lin_list, Lineage_v, mutaions_num, feature_mutations = feature_mut(
        download_path, lineage_file)
    print("The number of candidate parental lineages (Last updated date: Feb 08, 2023): ", len(lin_list), "\n")

    print("\n", "----------------------- STEP3: CovRecomb analysis -------------------------------------", "\n")
    print("The number of inputted genomes for CovRecomb pipline: ", len(Self_snp), "\n")
    print("The least number of the sequential feature mutations: ", str(len_UXY), "\n")
    print("The number of computational core used: ", str(cor_num), "\n")

    # the least number of sequential feature mutation must included in feature mutation pattern
    must_inA = "X" * len_UXY
    must_inB = "Y" * len_UXY

    # the output file
    col_names = [
        'sample_id','lineage_X', 'lineage_Y', \
    'mutation_pattern', "raw_p_value","adjusted_p_value","X_mutations", "Y_mutations", "shared_mutations", "denovo_mutations"
    ]
    with open(output_file, "w") as file_epi:
        for c in col_names[0:-1]:
            file_epi.write(c + ",")
        file_epi.write(col_names[-1] + "\n")

    os.chdir(os.getcwd() + "/")

    if len(Self_snp) <= cor_num:
        recombination_detection(Self_snp, Self_variants, Lineage_v, lin_list, must_inA,must_inB,feature_mutations, output_file, mutaions_num,max_bk_num,len_UXY)
    else:
        for i in range(1, cor_num + 1):
            globals()["p"+str(i)] = Process(target=recombination_detection,\
                args=(Self_snp[(i-1)*int(len(Self_snp)/cor_num):i*int(len(Self_snp)/cor_num)],\
                    Self_variants, Lineage_v, lin_list, must_inA,must_inB,feature_mutations, output_file, mutaions_num,max_bk_num,len_UXY))
                
        for i in range(1, cor_num + 1):
            eval("p" + str(i)).start()

        for i in range(1, cor_num + 1):
            eval("p" + str(i)).join()

    df = pd.read_csv(output_file)
    print("\n", "----------------------- Analysis completed --------------------------------------------", "\n")
    print("The number of putative inter-lineage recombinants: ", str(df.shape[0]), "\n")
    print("The total elapsed time：%.8s s" % (time.time() - maintime), "\n")


if __name__ == '__main__':
    main()
