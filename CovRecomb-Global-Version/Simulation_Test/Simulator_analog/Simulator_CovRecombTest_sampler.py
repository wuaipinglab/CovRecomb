'''
Author: your name
Date: 2022-04-11 16:55:30
LastEditTime: 2022-07-07 22:07:07
LastEditors: Sonia-Ljy lijysunny@sina.com
Description: 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
FilePath: /undefined/home/soniali/Desktop/03_recom_0308/10_verification/lineage_zhenshi_diff/compare_3seq_sampler/seq3_write_sampler.py
'''

from numpy.random import default_rng
import csv
from numpy import *
import pandas as pd
import argparse
import sys
import os

os.chdir(os.getcwd() + "/")


def lin_initial_sequence(pop_size, length=1000, ref=False):

    def generate_sequence(length):
        """
        This method will generate a sequence, and set the Sequence object's 
        sequence to that sequence.
        """
        sequence = ''
        for i in range(length):
            import random
            letter = random.choice(['A', 'T', 'G', 'C'])
            sequence += letter
        return sequence

    sequence_linA = generate_sequence(length)
    lin_seq = {}
    if ref == False:
        num = 0
        for i in range(pop_size):
            num += 1
            lin_ID = "lin" + str(num)

            lin_seq[lin_ID] = sequence_linA
    elif ref == True:
        lin_seq["reference"] = sequence_linA

    return lin_seq


def creat_dir(turns_file):
    import os
    if not os.path.exists(turns_file):
        os.makedirs(turns_file)


def choose_new_letter(letter):
    """
    This function chooses a new letter from ATGC that is
    different from the letter passed into the function.

    INPUTS:
    -	CHAR: letter
            the letter that will not be chosen from ATGC.
    """
    possible_letters = set(['A', 'T', 'G', 'C'])
    import random
    new_letter = random.choice(list(possible_letters.difference(set(letter))))  # 随机选取除了原碱基之外的一个碱基

    return new_letter


def mutate(sequence, mutation_site, qmatrix):
    '''
        A function that mutates a selected individual sequence and checks
        if the mutation was synonymous or nonsynonymous. Individual fitness
        gets decreased for each nonsynonymous mutation.

        Parameters:
            individual (dict): A dictionary representing an individual in
                               a population.
            mutation_sites (list): Sites where mutations occur.
            positions (dict): Dictionary with information about each
                                  coding position in a genome.

        Returns:
            individual (dict): Returns individual dictionary with a
                               mutated sequence.
    '''
    rng = default_rng()
    site = int(mutation_site)
    mutated_base = sequence[site]
    probs = qmatrix[mutated_base]
    possibleBases = list(probs.keys())
    probBases = list(probs.values())
    new_base = rng.choice(possibleBases, p=probBases)

    # new_base = choose_new_letter(mutated_base)
    return new_base


# Initialize sequences
def generate_ref_and_lin(seq_length, lineage_num, lin_diff_num):
    ref_seq = lin_initial_sequence(1, seq_length, ref=True)
    all_pop = {}
    mut_perlin = lin_diff_num

    all_pos = [x for x in range(0, seq_length)]
    # generate the seed (lineages)
    already_pos = []
    for seed_num in range(1, lineage_num + 1):
        seed_name = "lin" + str(seed_num)
        left_pos = [x for x in all_pos if x not in already_pos]
        import random
        positions = random.sample(left_pos, mut_perlin)
        seed_seq = ""
        for position in all_pos:
            if position in positions:
                letter = choose_new_letter(ref_seq["reference"][position])
                seed_seq += letter
            else:
                letter = ref_seq["reference"][position]
                seed_seq += letter

        all_pop[seed_name] = seed_seq
        # 选出已经发生过突变的位点
        already_pos += positions

    return ref_seq, all_pop


# 开始测试
def output_fv_file(generations, lin_fv75, turns_file):
    with open(turns_file + "/fv_75.txt", "w") as f:
        for l in lin_fv75:
            f.write(l + "," + str(generations) + ",")
            for v in lin_fv75[l]:
                f.write(v + ",")
            f.write("\n")


def extract_fv(DIRPATH):
    lineage_file = 'fv_75_norm_cluster75.txt'
    Lineage_v = {}
    linA_list = []
    feature_mutations = []
    with open(DIRPATH + lineage_file, 'r') as f:
        for i in csv.reader(f):
            if len(i[2:]) >= 1:
                linA_list.append(i[0])
                Lineage_v[i[0]] = i[2:]
                for v in i[2:]:
                    if v not in feature_mutations:
                        feature_mutations.append(v)
    return Lineage_v, linA_list, feature_mutations


def CovRecom_test(sample_child_mutate, len_UAB, Lineage_v, linA_list, feature_mutations, mutaions_num, code_path):
    import os
    os.chdir(code_path)
    import CovRecomb_method

    Strain_list_child = []
    variants_all_child = {}
    for ch in sample_child_mutate:
        Strain_list_child.append(ch)
        variants_all_child[ch] = sample_child_mutate[ch]

    child_num, child_correct_AB_num, child_false_AB_num, child_correct_AA_num, child_false_AA_num, child_zhenyin, child_jiayang, child_jiayin, child_zhenyang\
        = CovRecomb_method.recombination_detection_test(Strain_list_child, len_UAB, linA_list, variants_all_child, feature_mutations, Lineage_v, mutaions_num)

    return child_num, child_correct_AB_num, child_false_AB_num, child_correct_AA_num, child_false_AA_num, child_zhenyin, child_jiayang, child_jiayin, child_zhenyang


def extract_ref_mut(ref_seq, seq_samples, lin_all_mutate, s):
    for g in range(len(seq_samples[s])):
        if seq_samples[s][g] != ref_seq["reference"][g]:
            mut = str(g + 1) + "_" + str(seq_samples[s][g])
            lin_all_mutate.append(mut)
    return lin_all_mutate


def get_lin_name(x1_lin):

    if "bk" in x1_lin:
        lin_name = ""
        flag = 0
    elif x1_lin.count("_") == 0:
        lin_name = x1_lin
        flag = 1
    else:
        lin_name = "lin" + x1_lin.split("lin")[1]
        flag = 1

    return lin_name, flag


def generate_recom_bk(seq_length, generation, x1_lin, x2_lin, pop):
    import random
    breakup = random.sample(range(0, seq_length), 1)  # The random sampled breakpoint
    new_recom_name = "gen" + str(generation) + "_" + x1_lin + "bk" + str(breakup[0]) + x2_lin + "/"
    # add seq:
    left = pop[x1_lin][:breakup[0]]  # not including break point
    right = pop[x2_lin][breakup[0]:]  # including break point
    recom_seq = left + right
    pop[new_recom_name] = recom_seq  # new recombinant
    return pop


def test_ld_genr(num_seeds, lin_diff_num_list, turns_list, generations_list, temp_file, code_path, qmatrix, sample_rate, len_UAB=4, seq_length=29903, mut_rate=float(8E-4 / 52), rec_rate=0.1):
    print('process_id: (%s)...' % (os.getpid()))
    rng = default_rng()
    result_record = {}

    linA_list = []
    for i in range(1, num_seeds + 1):
        linA_list.append("lin" + str(i))

    for lin_diff_num in lin_diff_num_list:
        run_total_number = 0
        df_cor_rate_record = pd.DataFrame(columns=["lin_diff_num", "genera", "turn_times", "mean_AB_cr", "mean_AA_cr", "mean_FPR", "mean_TPR"])

        for genera in generations_list:
            # genera_number = 0
            df_result_record = pd.DataFrame(columns=[
                "lin_diff_num", "genera", "turns", "child_num", "child_correct_AB_num", "child_false_AB_num", "all_recom_correct_rate", "child_correct_AA_num", "child_false_AA_num", "all_yin_correct_rate", "child_zhenyin_d", "child_jiayang_b", "child_jiayin_c",
                "child_zhenyang_a", "x_FPR", "y_TPR", "usetime"
            ])

            for turns in range(0, turns_list):
                # turns = 0
                # Start for the turns number of independent trail
                print("lin_diff_num:", lin_diff_num, "  generation:", genera, "turns:", turns, "sample_rate:", sample_rate)
                # generate the reference seq and the targed number of lineage's sequence
                ref_seq, all_pop = generate_ref_and_lin(seq_length, num_seeds, lin_diff_num)

                # Output file name
                tagRun = "ld" + str(lin_diff_num) + "gene" + str(genera) + "turns" + str(turns)

                seq_samples = []
                import copy
                pop = copy.deepcopy(all_pop)

                # Loop for every generation
                for generation in range(1, genera + 1):
                    # Mutation
                    pool = copy.deepcopy(pop)
                    for lin in pool:
                        if "bk" not in lin:
                            new_lin_name = "gen" + str(generation) + "_" + lin

                            seed_sequence = pop[lin]
                            # The mutated sequence based on the seed (or former) sequence
                            mutations_seq = ""

                            # Calculated the number of mutate sites and positions based on the length of genomes.
                            poisson = rng.poisson(mut_rate, seq_length)  # rate is per site
                            num_mut = round(sum(poisson))

                            if num_mut > 0:
                                import random
                                positions = random.sample(range(0, seq_length), num_mut)  # 随机采出的突变位点

                                for position in range(0, len(seed_sequence)):
                                    if position in positions:
                                        letter = mutate(seed_sequence, position, qmatrix)
                                        mutations_seq += letter
                                    else:
                                        letter = seed_sequence[position]
                                        mutations_seq += letter
                                pop[new_lin_name] = mutations_seq
                            else:
                                pop[new_lin_name] = seed_sequence

                    # Recombination
                    poisson = rng.poisson(rec_rate, len(pop))  # len(pop)) #泊松分布 # rate is per individual
                    numRec = round(sum(poisson))
                    if numRec == 0:  # Enable at least one inter-lineage recombinant and one intra-lineage recombination in each generation
                        numRec = 1

                    already_id = list(pop.keys())

                    # Inter-lineage recombination
                    for i in range(numRec):
                        flag_x1 = flag_x2 = 0
                        lin_name_x1 = lin_name_x2 = ""
                        # Select the sequences sourcing from different seeds(lineages), should exclude recombinant.
                        while flag_x1 == 0 or flag_x2 == 0 or lin_name_x1 == lin_name_x2:
                            import random
                            x1_lin, x2_lin = random.sample(already_id, 2)
                            lin_name_x1, flag_x1 = get_lin_name(x1_lin)
                            lin_name_x2, flag_x2 = get_lin_name(x2_lin)

                        import random
                        pop = generate_recom_bk(seq_length, generation, x1_lin, x2_lin, pop)

                    # Intra-lineage recombination
                    for i in range(numRec):
                        flag_xAA = 0
                        while flag_xAA == 0:
                            xAA1_lin = random.sample(already_id, 1)[0]
                            lin_name_xAA, flag_xAA = get_lin_name(xAA1_lin)
                        # find a sequence souced from the same seed(lineage)
                        random.shuffle(already_id)
                        for aa in already_id:
                            if ("bk" not in aa) and (lin_name_xAA in aa) and (aa != xAA1_lin):
                                xAA2_lin = aa
                                break
                        pop = generate_recom_bk(seq_length, generation, xAA1_lin, xAA2_lin, pop)

                # The end of the simulation process and the simulted dataset in this independent trial has been generated.
                # 构建虚拟序列结束
                sampler_num = int(len(pop) / sample_rate)
                seq_samples_id = random.sample(pop.keys(), sampler_num)
                seq_samples = {}
                for key in pop:
                    if key in seq_samples_id:
                        seq_samples[key] = pop[key]

                # Extract feature mutations for each lineage (use 0.75 as cutoff value)
                lin_fv75 = {}
                for num in range(1, num_seeds + 1):
                    lin_all_mutate = []  # Record all the mutations within each lineage
                    # Record the number of samples within each lineage, except recombinants
                    count_target_lin = 0

                    for s in seq_samples:
                        if "bk" not in s:
                            if int(s.split("lin")[1]) == num:
                                count_target_lin += 1
                                # feature mutations
                                lin_all_mutate = extract_ref_mut(ref_seq, seq_samples, lin_all_mutate, s)
                                # parental lineages
                                seq_file0 = temp_file + "seq_method/"
                                creat_dir(seq_file0)
                                seq_file1 = seq_file0 + "lin_diff_num" + str(lin_diff_num)
                                creat_dir(seq_file1)
                                seq_file2 = seq_file1 + "/turns" + str(turns)
                                creat_dir(seq_file2)
                                with open(seq_file2 + "/lin_diff_num" + str(lin_diff_num) + "gener" + str(genera) + "turns" + str(turns) + "_parental_seq.fasta", "a+") as h:
                                    h.write(">" + s + "\n")
                                    h.write(seq_samples[s] + "\n")

                    fv75 = []
                    for m in set(lin_all_mutate):
                        if lin_all_mutate.count(m) >= (0.75 * count_target_lin):
                            fv75.append(m)

                    lin_fv75["lin" + str(num)] = fv75

                # Extract feature mutations for each recombinant
                sample_child_mutate = {}
                for r in seq_samples:
                    if "bk" in r:
                        temp_recom_mut = []
                        temp_recom_mut = extract_ref_mut(ref_seq, seq_samples, temp_recom_mut, r)
                        sample_child_mutate[r] = temp_recom_mut  # Mutations for each recombinant
                        # Output the fasta file for recombinants
                        with open(seq_file2 + "/lin_diff_num" + str(lin_diff_num) + "gener" + str(genera) + "turns" + str(turns) + "_recom_seq.fasta", "a+") as h:
                            h.write(">" + r + "\n")
                            h.write(seq_samples[r] + "\n")

                # part1  CovRecomb to detect recombinants
                feature_mutation = []
                for l in linA_list:
                    for v in lin_fv75[l]:
                        feature_mutation.append(v)
                feature_mutations = list(set(feature_mutation))
                mutaions_num = len(feature_mutations)

                import datetime
                starttime = datetime.datetime.now()
                child_num, child_correct_AB_num, child_false_AB_num, child_correct_AA_num, child_false_AA_num, child_zhenyin, child_jiayang, child_jiayin, child_zhenyang\
                    = CovRecom_test(sample_child_mutate, len_UAB, lin_fv75, linA_list, feature_mutations, mutaions_num, code_path)
                endtime = datetime.datetime.now()
                usetime = (endtime - starttime).microseconds

                if (child_correct_AB_num == child_false_AB_num == 0) or (child_correct_AA_num == child_false_AA_num == 0):
                    turns_list += 1
                    continue
                else:
                    # genera_number+=1
                    correct_rate = round(child_correct_AB_num / (child_correct_AB_num + child_false_AB_num), 4)
                    all_yin_correct_rate = round(child_correct_AA_num / (child_correct_AA_num + child_false_AA_num), 4)  # ZeroDivisionError: division by zero
                    if child_jiayang == child_zhenyin == 0:
                        x_FPR = "na"
                    else:
                        x_FPR = round(child_jiayang / (child_jiayang + child_zhenyin), 4)

                    y_TPR = round(child_zhenyang / (child_zhenyang + child_jiayin), 4)

                    result_record[tagRun] = [
                        lin_diff_num, genera, turns, child_num, child_correct_AB_num, child_false_AB_num, correct_rate, child_correct_AA_num, child_false_AA_num, all_yin_correct_rate, child_zhenyin, child_jiayang, child_jiayin, child_zhenyang, x_FPR, y_TPR, usetime
                    ]

                    df_result_record.loc[tagRun] = result_record[tagRun]

            # End of all the independent trails (turns)
            recom_file = temp_file + "covrecom/"
            creat_dir(recom_file)
            df_result_record.to_csv(recom_file + "ld" + str(lin_diff_num) + "gener" + str(genera) + "_truns" + str(turns + 1) + "_result" + ".csv")

            run_total_number += 1
            rr_geneer_yang_cr = df_result_record["all_recom_correct_rate"].tolist()
            rr_geneer_yin_cr = df_result_record["all_yin_correct_rate"].tolist()
            x_FPR_turns_mean = df_result_record["x_FPR"].tolist()
            y_TPR_turns_mean = df_result_record["y_TPR"].tolist()

            df_cor_rate_record.loc[run_total_number] = [lin_diff_num, genera, str(turns + 1), mean(rr_geneer_yang_cr), mean(rr_geneer_yin_cr), mean(x_FPR_turns_mean), mean(y_TPR_turns_mean)]

        df_cor_rate_record.to_csv(recom_file + "/lin_diff_num" + str(lin_diff_num) + "_result" + ".csv")


def main(sysargs=sys.argv[1:]):
    parser = argparse.ArgumentParser(description='Simulator_CovRecombTest', usage='''python3 Simulator_CovRecombTest_sampler.py -ld [15,21,38,46] -gen [8]''')

    parser.add_argument("-ld", help="The list number of differential mutations between simulated lineages, such as [46]", default=[46])
    parser.add_argument("-gen", "--generation", help="The list generation number for each simulation dataset \n Default: [4]", default=[8])
    parser.add_argument("-sd", "--seed_number", help="The number of seed (initial lineages) \n Default: 50", type=int, default=50)
    parser.add_argument("-sr", "--sample_rate", help="How many number of recombinants will be sampled for one sequence ", type=int, default=100)

    parser.add_argument("-turns", "--turns_number", help="The number of independent simulation trials \n Default: 10", default=10)
    parser.add_argument("-threshold", "--least_number_of_sequential_feature_mutations", help="The threthold for the least_number_of_sequential_feature_mutation \n Default: 4", default=4)
    parser.add_argument("-seqlength", "--seq_length", help="The length of ATCG for each simulated sequence \n Default: 29903", default=29903)
    parser.add_argument("-mt", "--mutation_rate", help="The nucleotide mutation rate for each generation \n Default: float(8E-4/52)", default=float(8E-4 / 52))
    parser.add_argument("-rt", "--recombination_rate", help="The recombination rate for each generation \n Default: 0.01", default=0.01)
    parser.add_argument("-cp", "--code_path", default=os.getcwd() + "/", help="The code file")

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    lin_diff_num_list = eval(args.ld)
    generations_list = eval(args.generation)
    num_seeds = int(args.seed_number)
    sample_rate = int(args.sample_rate)
    code_path = args.code_path
    temp_file = code_path + "compare_3seq/"
    turns_list = int(args.turns_number)
    len_UAB = int(args.least_number_of_sequential_feature_mutations)
    seq_length = int(args.seq_length)
    mut_rate = float(args.mutation_rate)
    rec_rate = float(args.recombination_rate)

    print("lin_diff_num_list", lin_diff_num_list)
    print("generations_list:", generations_list)
    print("num_seeds:", num_seeds)
    print("sample_rate:", sample_rate)
    print("code_path:", code_path)
    print("rec_rate:", rec_rate)

    creat_dir(temp_file)
    cor_num = len(lin_diff_num_list)
    # Based on empirical n=1180 mutations counts from 7 countries, It is sourced from the record of a previous research (Saymon Akther, 2021).
    qmatrix = {'A': {"C": 0.1083, "G": 0.7000, "T": 0.1917}, 'C': {"A": 0.0475, "G": 0.0033, "T": 0.9492}, 'G': {"A": 0.2102, "C": 0.0931, "T": 0.6967}, 'T': {"A": 0.1025, "C": 0.795, "G": 0.1025}}

    print("\n", "----------------------- Parameters loaded -----------------------", "\n")
    from multiprocessing import Process
    for i in range(1, cor_num + 1):
        globals()["p" + str(i)] = Process(target=test_ld_genr, args=(num_seeds, [lin_diff_num_list[i - 1]], turns_list, generations_list, temp_file, code_path, qmatrix, sample_rate, len_UAB, seq_length, mut_rate, rec_rate))

    for i in range(1, cor_num + 1):
        eval("p" + str(i)).start()

    for i in range(1, cor_num + 1):
        eval("p" + str(i)).join()

    print("\n", '----------------------- The end of simulator and the CovRecomb test----------------------- ', "\n")


if __name__ == '__main__':
    main()
