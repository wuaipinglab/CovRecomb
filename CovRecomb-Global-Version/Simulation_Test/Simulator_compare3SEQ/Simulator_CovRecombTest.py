'''
Input: The empirical substitution rates for "ATCG" bases.

Parameters: (1) "-smp", The number of differential feature mutations in each seed's generation process; Default: 6.
            (2) "-pp", "--parallel_prop", The probability for sequences with homologous mutation(s); Default: 0.1.
            (3) "-sgh", "--seed_gene_homo", The number of homologous mutation(s) in each seed generation process; Default: 1.
            (4) "-sr", "--sample_rate", The sample rate for all the generated sequences for analyzing; Default: 1.
            (5) "-gen", "--generation", The generation number for each simulation dataset in the format of a list; Default: [4,6].
            (6) "-sg", "--seed_gen", The generation number in lineage(seed) generation process; Default: 5.
            (7) "-sd", "--seed_number", The number of seed (initial lineages); Default: 5.
            (8) "-turns", "--turns_number", The number of independent simulation trials; Default: 10.
            (9) "-threshold", "--least_number_of_sequential_feature_mutations", The threthold for the least_number_of_sequential_feature_mutation; Default: 4.
            (10) "-seqlength", "--seq_length", The length of ATCG for each simulated sequence; Default: 29903", default=29903.
            (11) "-mt", "--mutation_rate", The nucleotide mutation rate for each generation; Default: float(8E-4/6).
            (12) "-rt", "--recombination_rate", The recombination rate for each generation; Default: 0.1.
            (13) "-comp", "--compare3SEQ", Whether to write down the fasta files for 3SEQ analyzing; Default: True.
            (14) "-f", "--temp_file", The file address; Default=os.getcwd()+"/".
            
Function: By simulating the viral transmission, mutation, and recombination process, all the generated the simulation sequences were extracted to test the coverage rate of CovRecomb.

Output: (1) Two .fasta files respectively representing the simulated non-recombinant sequences (_parental_seq.fasta) or the recombinant genomes (_recom_seq.fasta).
		(2) The detected coverage rate and the elapsed time of each trial for CovRecomb to detect inter-lineage recombinants (In "covrecom" folder).
'''

import random
from numpy.random import default_rng
import csv
import os
from numpy import *
import pandas as pd
from Bio import SeqIO
import random
import argparse
import sys
# os.chdir(os.getcwd()+"/")

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


def lin_initial_sequence(pop_size,temp_file,length = 1000,ref = True):
    with open(temp_file+"EPI_ISL_402125.fasta", "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            sequence_linA = str(record.seq)
            
    lin_seq = {}
    if ref == False:
        num = 0
        for i in range(pop_size):
            num += 1
            lin_ID = "lin"+ str(num)
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
    new_letter = random.choice(list(
        possible_letters.difference(set(letter))))
    return new_letter


rng = default_rng() # Random number generator
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
                               mutated sequence and changed fitness
                               value for each nonsynonymous mutation.
    '''
    site = int(mutation_site)
    mutated_base = sequence[site]
    probs = qmatrix[mutated_base]
    possibleBases = list(probs.keys())
    probBases = list(probs.values())
    new_base =  rng.choice(possibleBases, p = probBases)
    # new_base = choose_new_letter(mutated_base)
    return new_base


# Initialize sequences
def generate_ref_and_lin(seq_length, temp_file,parallel_prop, seed_gen, seed_mut_perlin,seed_gene_homo,num_seeds):
    ref_seq = lin_initial_sequence(1, temp_file,seq_length, ref = True)
    temp_pop = {}
    all_pos = [x for x in range(0, seq_length)]

    import random
    for gen in range(seed_gen):
        gen += 1
        homo_mut = random.sample(all_pos, seed_gene_homo)
        seq_homo = ""
        for site in all_pos:
            if site in homo_mut:
                seq_homo += choose_new_letter(ref_seq["reference"][site])
            else:
                seq_homo += ref_seq["reference"][site]

        mut_lin_num = int((2**gen)*parallel_prop) # The number of sequences with homologous mutations
        import random
        select_seq_with_homo = random.sample(range(2**gen), mut_lin_num)
        already_pos = []

        for n in range(2**gen):
            seed_name = "temp_gen"+str(gen)+"_seq"+str(n)
            left_pos = [x for x in all_pos if x not in already_pos]
            import random
            positions = random.sample(left_pos, seed_mut_perlin)
            already_pos.extend(positions)
            seed_seq = ""
            if gen >= 2:
                selected_seed_name = temp_pop["temp_gen"+str(int(gen-1))+"_seq"+str(int(n/2))]
            elif gen == 1:
                selected_seed_name = ref_seq["reference"]
                
            for position in all_pos:
                if position in positions:
                    letter = choose_new_letter(selected_seed_name[position])
                    seed_seq += letter
                elif position in homo_mut and n in select_seq_with_homo:
                    seed_seq += seq_homo[position] # homologous mutations
                else:
                    letter = selected_seed_name[position]
                    seed_seq += letter
            temp_pop[seed_name] = seed_seq

    last_pop = []
    for i in temp_pop:
        if "gen"+str(seed_gen) in i :
            last_pop.append(i)
    import random
    seed_name = random.sample(last_pop, num_seeds)
    seed_pop = {}
    for i in range(1,num_seeds+1):
        seed_pop["lin"+str(i)] = temp_pop[seed_name[i-1]]
        
    seed_diff = []
    for p in range(1,num_seeds+1):
        count = 0 
        for i in range(1,num_seeds+1):
            if p != i:
                for s in range(29903):
                    if seed_pop["lin"+str(p)][s] != seed_pop["lin"+str(i)][s]:
                        count+=1
                seed_diff.append(count)
    lin_diff_num = mean(seed_diff)

    linA_list = list(seed_pop.keys())
    return ref_seq, seed_diff, lin_diff_num, linA_list, seed_pop


# based on empirical n=1180 mutations counts from 7 countries
qmatrix = {
    'A': { "C": 0.1083, "G": 0.7000, "T": 0.1917 },
    'C': { "A": 0.0475, "G": 0.0033, "T": 0.9492 },
    'G': { "A": 0.2102, "C": 0.0931, "T": 0.6967 },
    'T': { "A": 0.1025, "C": 0.795, "G": 0.1025 }
    }

######开始测试
def output_fv_file(generations, lin_fv75, turns_file):
    with open(turns_file+"/fv_75.txt", "w") as f:
        for l in lin_fv75:
            f.write(l+","+str(generations)+",")
            for v in lin_fv75[l]:
                f.write(v+",")
            f.write("\n")


def extract_fv(DIRPATH):
    lineage_file = 'fv_75_norm_cluster75.txt'
    Lineage_v = {}
    linA_list = []
    feature_mutations = []
    with open(DIRPATH+lineage_file,'r') as f:
        for i in csv.reader(f):
            if len(i[2:])>= 1:
                linA_list.append(i[0])
                Lineage_v[i[0]] = i[2:]
                for v in i[2:]:
                    if v not in feature_mutations:
                        feature_mutations.append(v)
    return Lineage_v,linA_list,feature_mutations

def CovRecom_test(sample_child_mutate, len_UAB, Lineage_v, linA_list, feature_mutations, mutaions_num,temp_file):
    os.chdir(temp_file)
    import CovRecomb_method
    Strain_list_child = []
    variants_all_child  = {}
    for ch in sample_child_mutate:
        Strain_list_child.append(ch)
        variants_all_child[ch] = sample_child_mutate[ch]        
           
    # child_num,child_correct_AB_num,child_false_AB_num \
    child_num,child_correct_AB_num,child_false_AB_num,child_correct_AA_num,child_false_AA_num,child_zhenyin,child_jiayang,child_jiayin,child_zhenyang\
        = CovRecomb_method.recombination_detection_test(Strain_list_child,len_UAB,linA_list,variants_all_child,feature_mutations,Lineage_v,mutaions_num)             
    return child_num,child_correct_AB_num,child_false_AB_num,child_correct_AA_num,child_false_AA_num,child_zhenyin,child_jiayang,child_jiayin,child_zhenyang


def extract_ref_mut(ref_seq, seq_samples, lin_all_mutate, s):
    for g in range(len(seq_samples[s])): 
        if seq_samples[s][g] != ref_seq["reference"][g]:
            mut = str(g+1)+"_"+str(seq_samples[s][g])
            lin_all_mutate.append(mut)
    return lin_all_mutate


def get_lin_name(x1_lin):
    if "bk" in x1_lin:
        lin_name= ""
        flag = 0
    elif x1_lin.count("_") == 0:
        lin_name = x1_lin
        flag = 1
    else:
        lin_name = "lin"+x1_lin.split("lin")[1]
        flag = 1

    return lin_name,flag

def generate_recom_bk(seq_length,generation,x1_lin,x2_lin,pop):
    breakup = random.sample(range(0, seq_length), 1) # Random breakpoint
    new_recom_name = "gen"+str(generation)+"_"+x1_lin + "bk"+str(breakup[0])+x2_lin+"/"
    # add seq:
    left = pop[x1_lin][:breakup[0]] # not including break point
    right = pop[x2_lin][breakup[0]:] # including break point
    recom_seq = left + right
    pop[new_recom_name] = recom_seq # new recombinant
    return pop


qmatrix = {'A': { "C": 0.1083, "G": 0.7000, "T": 0.1917 },
        'C': { "A": 0.0475, "G": 0.0033, "T": 0.9492 },
        'G': { "A": 0.2102, "C": 0.0931, "T": 0.6967 },
        'T': { "A": 0.1025, "C": 0.795, "G": 0.1025 }}
    
column_names = ["num_seeds","generations_list","sampled_numer","turns_list","sample_rate","seed_gene_homo_rate","seed_mut_perlin","UAB","lin_diff_num","mean_lin_diff","FPR","FDR","TPR"]
with open("/home/soniali/Desktop/03_CovRecomb_review0819/simulation/simulation_renew/ROC_results.txt","a+") as f:
            f.write(",".join(column_names)+"\n")

df_cor_rate_record = pd.DataFrame(columns = ["lin_diff_num","genera","turn_times","mean_AB_cr","mean_AA_cr","mean_FPR","mean_FDR_wufaxianlv","mean_TPR","mean_usetime"] )

def simulator_withhomo(turns_list,seed_mut_perlin,generations_list, rec_rate, compare3SEQ,seq_length,temp_file,parallel_prop,seed_gen,seed_gene_homo,num_seeds,len_UAB,mut_rate,sample_rate):
    run_total_number = 0
    result_record = {}
    for genera in generations_list:
        import copy
        turns_list_deep = copy.deepcopy(turns_list)
        df_result_record = pd.DataFrame(columns = ["lin_diff_num","genera","turns","child_num","child_correct_AB_num","child_false_AB_num","all_recom_correct_rate",\
        "child_correct_AA_num","child_false_AA_num","all_yin_correct_rate",\
            "child_zhenyin_d","child_jiayang_b","child_jiayin_c","child_zhenyang_a",\
                "x_FPR","x_FDR_wufaxianlv","y_TPR","usetime"])
        lin_diff_num_list = []
        for turns in range(1,turns_list_deep+1):
            ## generate ref_seq and lineage_seq
            ref_seq, seed_diff, seed_diff_num, linA_list, all_pop = generate_ref_and_lin(seq_length, temp_file,parallel_prop, seed_gen, seed_mut_perlin,seed_gene_homo,num_seeds)
            ## start the No. test
            seq_samples = []
            import copy
            pop = copy.deepcopy(all_pop)
            # Loop for every process
            for generation in range(1, genera + 1):
                all_pos = [x for x in range(0,seq_length)]
                import random
                homo_mut = random.sample(all_pos, seed_gene_homo)
                seq_homo = ""
                for site in all_pos:
                    if site in homo_mut:
                        seq_homo += mutate(ref_seq["reference"],site,qmatrix)
                    else:
                        seq_homo += ref_seq["reference"][site]

                mut_lin_num = int(len(pop)*parallel_prop) 
                select_seq_with_homo = random.sample(range(len(pop)),mut_lin_num) 
                
                ## Mutation
                pool = copy.deepcopy(pop)
                count = 0
                for lin in pool:
                    count+=1
                    if "bk" not in lin:
                        new_lin_name = "gen"+str(generation)+"_"+lin
                        seed_sequence = pop[lin]
                        mutations_seq = "" # The mutated sequence based on the seed (or former) sequence
                        # Calculated the number of mutate sites and positions based on the length of genomes.
                        num_mut = round(sum(rng.poisson(mut_rate,seq_length))) # rate is per site

                        if num_mut > 0:
                            import random
                            positions = random.sample(range(0, seq_length), num_mut) 

                            if (count-1) not in select_seq_with_homo:
                                for position in range(0,len(seed_sequence)):
                                    if position in positions:
                                        letter = mutate(seed_sequence, position,qmatrix)
                                        mutations_seq += letter
                                    else:
                                        letter = seed_sequence[position] 
                                        mutations_seq += letter
                                        
                                pop[new_lin_name] = mutations_seq
                            else:
                                for position in range(0,len(seed_sequence)):
                                    if position in positions:
                                        letter = mutate(seed_sequence, position,qmatrix)
                                        mutations_seq += letter
                                    elif position in homo_mut:
                                        letter = seq_homo[position] 
                                        mutations_seq += letter
                                    else:
                                        letter = seed_sequence[position] 
                                        mutations_seq += letter
                        elif num_mut == 0:
                            for position in all_pos:
                                if position in homo_mut:
                                    letter = seq_homo[position]
                                    mutations_seq += letter
                                else:
                                    letter = seed_sequence[position] 
                                    mutations_seq += letter
                        pop[new_lin_name] = seed_sequence

                ## Recombination
                numRec = round(sum(rng.poisson(rec_rate, len(pop))))  
                if numRec == 0: # Enable at least one inter-lineage recombinant and one intra-lineage recombination in each generation
                    numRec = 1
                already_id = list(pop.keys())
                # Inter-lineage recombination 
                for i in range(numRec):
                    flag_x1 = flag_x2 = 0
                    lin_name_x1 = lin_name_x2 = ""
                    # Select the sequences sourcing from different seeds(lineages), should exclude recombinant.
                    while flag_x1 == 0 or  flag_x2 == 0 or lin_name_x1 == lin_name_x2:
                        import random
                        x1_lin, x2_lin = random.sample(already_id, 2)
                        lin_name_x1,flag_x1 = get_lin_name(x1_lin)
                        lin_name_x2,flag_x2 = get_lin_name(x2_lin)

                    import random
                    pop = generate_recom_bk(seq_length, generation, x1_lin, x2_lin, pop)

                # Intra-lineage recombination 
                for i in range(numRec):
                    flag_xAA = 0
                    while flag_xAA == 0:
                        xAA1_lin = random.sample(already_id, 1)[0]
                        lin_name_xAA,flag_xAA = get_lin_name(xAA1_lin)
                    # find a sequence souced from the same seed(lineage)
                    random.shuffle(already_id)
                    for aa in already_id:
                        if ("bk" not in aa) and (lin_name_xAA in aa) and (aa != xAA1_lin):
                            xAA2_lin = aa
                            break
                    pop = generate_recom_bk(seq_length,generation,xAA1_lin,xAA2_lin,pop)
            ## The end of the simulation process and the simulted dataset in this independent trial has been generated.
            sampler_num = int(len(pop)/sample_rate)
            seq_samples_id = random.sample(pop.keys(), sampler_num)
            seq_samples = {}
            for key in pop:
                if key in seq_samples_id:
                    seq_samples[key] = pop[key]
            
            sampled_numer = str(len(seq_samples))
            ## Extract feature mutations for each lineage (use 0.75 as cutoff value)
            lin_fv75 = {}
            usetime = 0
            for num in range(1,num_seeds+1):
                lin_all_mutate = [] # Record all the mutations within each lineage
                count_target_lin = 0 # Record the number of samples within each lineage, except recombinants
                for s in seq_samples:
                    if "bk" not in s and int(s.split("lin")[1]) == num:
                        count_target_lin += 1
                        # Get the feature mutations           
                        lin_all_mutate =  extract_ref_mut(ref_seq, seq_samples, lin_all_mutate, s)
                import datetime 
                starttime = datetime.datetime.now()
                fv75 = []
                for m in set(lin_all_mutate):
                    if lin_all_mutate.count(m) >= (0.75*count_target_lin): 
                        fv75.append(m)
                
                lin_fv75["lin"+str(num)] = fv75
                
                endtime = datetime.datetime.now()
                usetime = usetime + (endtime - starttime).microseconds
            
            lin_diff_num = []
            for lin in lin_fv75:
                for lin2 in lin_fv75:
                    lin_diff_num.append(len(set(lin_fv75[lin]) ^ set(lin_fv75[lin2])))
                    
            lin_diff_num = int(mean(lin_diff_num))
            lin_diff_num_list.append(lin_diff_num)
            ## the filename of the output file
            print("seed_diff_num:",str(seed_diff_num),"  lin_diff_num:", lin_diff_num,"  generation:",genera,"turns:",turns)
            tagRun ="s"+str(num_seeds)+"smpl"+str(seed_mut_perlin)+"gene"+str(genera)+"turns"+str(turns)+"RT"+str(rec_rate) # "len"+str(seq_length)+"sd"+str(num_seeds)+
            if compare3SEQ == True: 
                for s in seq_samples:
                    if "bk" not in s :
                        seq_file0 = temp_file+"seq_method/"
                        creat_dir(seq_file0)
                        seq_file1 = seq_file0+"gener"+str(genera)
                        creat_dir(seq_file1)
                        seq_file2 = seq_file1+"/turns"+str(turns)
                        creat_dir(seq_file2)
                        with open(seq_file2+"/gener"+str(genera)+"turns"+str(turns)+"_parental_seq.fasta","a+") as h:
                            h.write(">"+s+"\n")
                            h.write(seq_samples[s]+"\n")
                            
            ## Extract feature mutations for each recombinant
            sample_all_mutate = {}
            for r in seq_samples:
                if compare3SEQ == True and "bk" in r: 
                    temp_all_mut = []
                    temp_all_mut = extract_ref_mut(ref_seq, seq_samples, temp_all_mut, r)
                    sample_all_mutate[r] = temp_all_mut  ## Mutations for each recombinant
                    
                    ###输出重组子代库
                    with open(seq_file2+"/genera"+str(genera)+"turns"+str(turns)+"_recom_seq.fasta","a+") as h:
                        h.write(">"+r+"\n")
                        h.write(seq_samples[r]+"\n")
                elif compare3SEQ == False:
                    temp_all_mut = []
                    temp_all_mut = extract_ref_mut(ref_seq, seq_samples, temp_all_mut, r)
                    sample_all_mutate[r] = temp_all_mut  ## Mutations for each recombinant

            ############## part1  CovRecomb to detect recombinants
            starttime = datetime.datetime.now()
            feature_mutation = []
            for l in linA_list:
                for v in lin_fv75[l]:
                    feature_mutation.append(v)
            feature_mutations = list(set(feature_mutation))
            mutaions_num = len(feature_mutations)
            endtime = datetime.datetime.now()
            usetime = usetime + (endtime - starttime).microseconds
            
            import datetime
            starttime = datetime.datetime.now()
            child_num,child_correct_AB_num,child_false_AB_num,child_correct_AA_num,child_false_AA_num,child_zhenyin,child_jiayang,child_jiayin,child_zhenyang\
                = CovRecom_test(sample_all_mutate, len_UAB, lin_fv75, linA_list, feature_mutations, mutaions_num, temp_file)
            endtime = datetime.datetime.now()
            usetime = usetime + (endtime - starttime).microseconds
            
            if (child_correct_AB_num == child_false_AB_num == 0) or (child_correct_AA_num == child_false_AA_num == 0) or ((child_jiayang+child_zhenyin) == 0) or ((child_jiayang+child_zhenyang) == 0):
                turns_list_deep +=1
                continue
            else:
                correct_rate = round(child_correct_AB_num/(child_correct_AB_num+child_false_AB_num),4)
                all_yin_correct_rate = round(child_correct_AA_num/(child_correct_AA_num+child_false_AA_num),4)
                if child_jiayang == child_zhenyin == 0:
                    x_FPR = "na"
                    x_FDR_wufaxianlv = "na"
                else:
                    x_FPR = round(child_jiayang/(child_jiayang+child_zhenyin),4)
                    x_FDR_wufaxianlv = round(child_jiayang/(child_jiayang+child_zhenyang),4)
                    
                y_TPR = round(child_zhenyang/(child_zhenyang+child_jiayin),4)
                result_record[tagRun] = [lin_diff_num,genera,turns,child_num,child_correct_AB_num,child_false_AB_num,correct_rate,\
                    child_correct_AA_num,child_false_AA_num,all_yin_correct_rate,child_zhenyin,child_jiayang,child_jiayin,child_zhenyang,\
                    x_FPR,x_FDR_wufaxianlv, y_TPR,usetime]

                df_result_record.loc[tagRun] = result_record[tagRun]
        
        ### End of all the independent trails (turns)
        recom_file = temp_file+"CovRecomb_compare/"
        creat_dir(recom_file)
        df_result_record.to_csv(recom_file+"seed_mut_perlin_"+str(seed_mut_perlin)+"_ld_UAB"+str(len_UAB)+"_"+str(lin_diff_num)+"gener"+str(genera)+"_truns"+str(turns)+"_result"+".csv") 
        run_total_number+=1
        rr_geneer_yang_cr = df_result_record["all_recom_correct_rate"].tolist()
        rr_geneer_yin_cr = df_result_record["all_yin_correct_rate"].tolist()
        x_FPR_turns_mean = df_result_record["x_FPR"].tolist()
        x_FDR_wufaxianlv_turns_mean = df_result_record["x_FDR_wufaxianlv"].tolist()
        y_TPR_turns_mean = df_result_record["y_TPR"].tolist()
        use_time_mean = df_result_record["usetime"].tolist()
        lin_diff_num_mean = mean(lin_diff_num_list)

        df_cor_rate_record.loc[run_total_number] = [lin_diff_num_mean,genera,str(turns),mean(rr_geneer_yang_cr),\
            mean(rr_geneer_yin_cr),mean(x_FPR_turns_mean),mean(x_FDR_wufaxianlv_turns_mean),mean(y_TPR_turns_mean),mean(use_time_mean)]

    df_cor_rate_record.to_csv(recom_file+tagRun+"_result"+".csv")
    return tagRun, y_TPR_turns_mean, df_cor_rate_record
   


def main(sysargs=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        description='Simulator_CovRecombTest',
        usage='''python Simulator_CovRecombTest.py -smp 6 -sr 1''')

    parser.add_argument("-smp", help="The number of differential feature mutations in each seed's generation process \n Default: 6", default=6)
    parser.add_argument("-pp", "--parallel_prop", help="The probability for sequences with homologous mutation(s) \n Default: 0.1", default= 0.1)
    parser.add_argument("-sgh", "--seed_gene_homo", help="The number of homologous mutation(s) in each seed generation process\n Default: 1", default=1)
    parser.add_argument("-sr", "--sample_rate", help="The sample rate for all the generated sequences for analyzing \n Default: 1", default=1)
    parser.add_argument("-gen", "--generation", help="The generation number for each simulation dataset in the format of a list \n Default: [4,6]", default=[4,6])
    parser.add_argument("-sg", "--seed_gen", help="The generation number in lineage(seed) generation process \n Default: 5", default=5)
    parser.add_argument("-sd", "--seed_number", help="The number of seed (initial lineages) \n Default: 5", type=int, default=5)
    parser.add_argument("-turns", "--turns_number", help="The number of independent simulation trials \n Default: 10", default=10)
    parser.add_argument("-threshold", "--least_number_of_sequential_feature_mutations", help="The threthold for the least_number_of_sequential_feature_mutation \n Default: 4", default=4)
    parser.add_argument("-seqlength", "--seq_length", help="The length of ATCG for each simulated sequence \n Default: 29903", default=29903)
    parser.add_argument("-mt", "--mutation_rate", help="The nucleotide mutation rate for each generation \n Default: float(8E-4/6)", default=float(8E-4 / 6))
    parser.add_argument("-rt", "--recombination_rate", help="The recombination rate for each generation \n Default: 0.1", default=0.1)
    parser.add_argument("-comp", "--compare3SEQ", help="Whether to write down the fasta files for 3SEQ analyzing \n Default: True", default=True)
    parser.add_argument("-f", "--temp_file", help="The file address \n", default=os.getcwd()+"/")
    
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    
    seed_mut_perlin = int(args.smp)
    seed_gen = eval(args.seed_gen)
    generations_list = eval(args.generation)
    num_seeds = int(args.seed_number)
    turns_list = int(args.turns_number)
    len_UAB = int(args.least_number_of_sequential_feature_mutations)
    seq_length = int(args.seq_length)
    mut_rate = float(args.mutation_rate)
    rec_rate = float(args.recombination_rate)
    seed_gene_homo = int(args.seed_gene_homo)
    sample_rate = int(args.sample_rate)
    parallel_prop = float(args.parallel_prop)
    compare3SEQ = args.compare3SEQ
    temp_file = args.temp_file
    creat_dir(temp_file)

    tagRun, y_TPR_turns_mean,df_cor_rate_record = simulator_withhomo(turns_list,seed_mut_perlin,\
        generations_list, rec_rate, compare3SEQ,seq_length,temp_file,parallel_prop,seed_gen,\
            seed_gene_homo,num_seeds,len_UAB,mut_rate,sample_rate)
    print(tagRun, y_TPR_turns_mean,df_cor_rate_record)

if __name__ == '__main__':
    main()
