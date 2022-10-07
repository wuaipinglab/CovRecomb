

from numpy.random import default_rng
import csv
from numpy import *
import pandas as pd
import os
os.chdir(os.getcwd()+"/")

def lin_initial_sequence(pop_size,length = 1000,ref = False):
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


def mutate(sequence, mutation_site,qmatrix):
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
    new_base =  rng.choice(possibleBases, p = probBases)
    # new_base = choose_new_letter(mutated_base)
    return new_base


# Initialize sequences
def generate_ref_and_lin(seq_length, parallel_prop, seed_gen, seed_mut_perlin,seed_gene_homo_rate):
    ref_seq = lin_initial_sequence(1, seq_length, ref = True)
    all_pop = {}
    temp_pop = {}
    all_pos = [x for x in range(0,seq_length)]

    import random
    for gen in range(seed_gen):
        gen += 1
        homo_mut = random.sample(all_pos, int(seed_mut_perlin*seed_gene_homo_rate))
        seq_homo = ""
        for site in all_pos:
            if site in homo_mut:
                seq_homo += choose_new_letter(ref_seq["reference"][site])
            else:
                seq_homo += ref_seq["reference"][site]

        mut_lin_num = int((2**gen)*parallel_prop)
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
                    seed_seq += seq_homo[position]
                else:
                    letter = selected_seed_name[position]
                    seed_seq += letter
            temp_pop[seed_name] = seed_seq

    seed_diff = []
    for p in range((2**seed_gen)):
        count = 0 
        for i in range(29903):
            if temp_pop["temp_gen"+str(seed_gen)+"_seq0"][i] != temp_pop["temp_gen"+str(seed_gen)+"_seq"+str(p)][i]:
                count+=1
        seed_diff.append(count)
        
    all_pop = {}
    num = 0
    for seed in temp_pop:
        if "gen"+str(seed_gen) in seed:
            num += 1
            all_pop["lin"+str(num)] = temp_pop[seed]
    
    linA_list = list(all_pop.keys())
    return ref_seq, seed_diff, linA_list, all_pop


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


def CovRecom_test(sample_all_mutate, len_UAB, Lineage_v, linA_list, feature_mutations, mutaions_num,code_path):
    # Lineage_v = lin_fv75
    import os
    os.chdir(code_path)
    import CovRecomb_method
    
    Strain_list_all = []
    variants_all_all = {}
    for ch in sample_all_mutate:
        Strain_list_all.append(ch)
        variants_all_all[ch] = sample_all_mutate[ch]        
           
    child_num,child_correct_AB_num,child_false_AB_num,child_correct_AA_num,child_false_AA_num,child_zhenyin,child_jiayang,child_jiayin,child_zhenyang\
        = CovRecomb_method.recombination_detection_test(Strain_list_all,len_UAB,linA_list,variants_all_all,feature_mutations,Lineage_v,mutaions_num)
    return child_num,child_correct_AB_num,child_false_AB_num,child_correct_AA_num,child_false_AA_num,child_zhenyin,child_jiayang,child_jiayin,child_zhenyang


def extract_ref_mut(ref_seq, seq_samples, lin_all_mutate, s):
    for g in range(len(seq_samples[s])):
        if seq_samples[s][g] != ref_seq["reference"][g]:
            mut = str(g+1)+"_"+str(seq_samples[s][g])
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
        lin_name = "lin"+x1_lin.split("lin")[1]
        flag = 1
    return lin_name, flag


def generate_recom_bk(seq_length,generation,x1_lin,x2_lin,pop):
    import random
    breakup = random.sample(range(0, seq_length), 1) # The random sampled breakpoint
    new_recom_name = "gen"+str(generation)+"_"+x1_lin + "bk"+str(breakup[0])+x2_lin+"/"
    # add seq:
    left = pop[x1_lin][:breakup[0]] # not including break point
    right = pop[x2_lin][breakup[0]:] # including break point
    recom_seq = left + right
    pop[new_recom_name] = recom_seq # new recombinant
    return pop


def test_ld_genr(homo_num,seed_gen, seed_gene_homo_rate, parallel_prop, seed_mut_perlin, turns_list, generations_list, temp_file, code_path, qmatrix, sample_rate, len_UAB = 4, seq_length = 29903, mut_rate = float(8E-4/52), rec_rate = 0.027):
    rng = default_rng()
    result_record = {}
    run_total_number = 0
    df_cor_rate_record = pd.DataFrame(columns = ["lin_diff_num","genera","turn_times","mean_AB_cr","mean_AA_cr","mean_FPR","mean_FDR_wufaxianlv","mean_TPR"] )
    df_result_record = pd.DataFrame(columns = ["lin_diff_num","genera","turns","child_num","child_correct_AB_num","child_false_AB_num","all_recom_correct_rate",\
            "child_correct_AA_num","child_false_AA_num","all_yin_correct_rate",\
                "child_zhenyin_d","child_jiayang_b","child_jiayin_c","child_zhenyang_a",\
                    "x_FPR","x_FDR_wufaxianlv","y_TPR"])
    
    for genera in generations_list:
        lin_diff_num_list = []
        for turns in range(1,turns_list+1):
            # Start for the turns number of independent trail
            print("seed_gen:",str(seed_gen),"seed_mut_perlin:",str(seed_mut_perlin),"UAB:",str(len_UAB),"turns:",turns)
            # generate the reference seq and the targed number of lineage's sequence
            ref_seq, seed_diff, linA_list, all_pop = generate_ref_and_lin(seq_length,parallel_prop, seed_gen, seed_mut_perlin,seed_gene_homo_rate)
            lin_diff_num = str(seed_diff[1])+"_"+str(max(seed_diff))
            # Output file name
            tagRun ="ld"+str(lin_diff_num)+"gene"+str(genera)+"turns"+str(turns)
            
            seq_samples = []
            import copy
            pop = copy.deepcopy(all_pop)

            # Loop for every generation
            for generation in range(1, genera + 1):
                all_pos = [x for x in range(0,seq_length)]
                import random
                homo_mut = random.sample(all_pos, homo_num)
                seq_homo = ""
                for site in all_pos:
                    if site in homo_mut:
                        seq_homo += choose_new_letter(ref_seq["reference"][site])
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
                        poisson = rng.poisson(mut_rate,seq_length) # rate is per site
                        num_mut = round(sum(poisson))

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
                poisson = rng.poisson(rec_rate, len(pop))
                numRec = round(sum(poisson))
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
            num_seeds = 2**seed_gen
            usetime = 0
            for num in range(1,num_seeds+1):
                lin_all_mutate = [] # Record all the mutations within each lineage
                count_target_lin = 0 # Record the number of samples within each lineage, except recombinants
                for s in seq_samples:
                    if "bk" not in s and int(s.split("lin")[1][0]) == num:
                        count_target_lin += 1
                        # Extract feature mutations                     
                        lin_all_mutate =  extract_ref_mut(ref_seq, seq_samples, lin_all_mutate, s)
                fv75 = []
                for m in set(lin_all_mutate):
                    if lin_all_mutate.count(m) >= (0.75*count_target_lin): 
                        fv75.append(m)
                
                lin_fv75["lin"+str(num)] = fv75

            
            lin_diff_num = []
            for lin in lin_fv75:
                for lin2 in lin_fv75:
                    lin_diff_num.append(len(set(lin_fv75[lin]) ^ set(lin_fv75[lin2])))
                    
            lin_diff_num = int(mean(lin_diff_num))
            lin_diff_num_list.append(lin_diff_num)
           
            ## Extract feature mutations for each recombinant
            sample_all_mutate = {}
            for r in seq_samples:
                temp_all_mut = []
                temp_all_mut = extract_ref_mut(ref_seq, seq_samples, temp_all_mut, r)
                sample_all_mutate[r] = temp_all_mut  ## Mutations for each recombinant

            ############## part1  CovRecomb to detect recombinants
            feature_mutation = []
            for l in linA_list:
                for v in lin_fv75[l]:
                    feature_mutation.append(v)
            feature_mutations = list(set(feature_mutation))
            mutaions_num = len(feature_mutations)

            child_num,child_correct_AB_num,child_false_AB_num,child_correct_AA_num,child_false_AA_num,child_zhenyin,child_jiayang,child_jiayin,child_zhenyang\
                = CovRecom_test(sample_all_mutate, len_UAB, lin_fv75, linA_list, feature_mutations, mutaions_num,code_path)
            
            if (child_correct_AB_num == child_false_AB_num == 0) or (child_correct_AA_num == child_false_AA_num == 0) or ((child_jiayang+child_zhenyin) == 0) or ((child_jiayang+child_zhenyang) == 0):
                turns_list +=1
                continue
            else:
                correct_rate = round(child_correct_AB_num/(child_correct_AB_num+child_false_AB_num),4)
                all_yin_correct_rate = round(child_correct_AA_num/(child_correct_AA_num+child_false_AA_num),4) ### ZeroDivisionError: division by zero
                if child_jiayang == child_zhenyin == 0:
                    x_FPR = "na"
                    x_FDR_wufaxianlv = "na"
                else:
                    x_FPR = round(child_jiayang/(child_jiayang+child_zhenyin),4)
                    x_FDR_wufaxianlv = round(child_jiayang/(child_jiayang+child_zhenyang),4)
                    
                y_TPR = round(child_zhenyang/(child_zhenyang+child_jiayin),4)
                result_record[tagRun] = [lin_diff_num,genera,turns,child_num,child_correct_AB_num,child_false_AB_num,correct_rate,\
                    child_correct_AA_num,child_false_AA_num,all_yin_correct_rate,child_zhenyin,child_jiayang,child_jiayin,child_zhenyang,\
                    x_FPR,x_FDR_wufaxianlv, y_TPR]

                df_result_record.loc[tagRun] = result_record[tagRun] #### 
        
        # End of all the independent trails (turns)
        recom_file = temp_file+"CovRecom/"
        creat_dir(recom_file)
        # df_result_record.to_csv(recom_file+"seed_mut_perlin_"+str(seed_mut_perlin)+"_ld_UAB"+str(len_UAB)+"_"+str(lin_diff_num)+"gener"+str(genera)+"_truns"+str(turns+1)+"_result"+".csv") 
        run_total_number+=1
        rr_geneer_yang_cr = df_result_record["all_recom_correct_rate"].tolist()
        rr_geneer_yin_cr = df_result_record["all_yin_correct_rate"].tolist()
        x_FPR_turns_mean = df_result_record["x_FPR"].tolist()
        x_FDR_wufaxianlv_turns_mean = df_result_record["x_FDR_wufaxianlv"].tolist()
        y_TPR_turns_mean = df_result_record["y_TPR"].tolist()
        mean_lin_diff = mean(lin_diff_num_list)

        df_cor_rate_record.loc[run_total_number] = [str(mean_lin_diff),genera,str(turns+1),mean(rr_geneer_yang_cr),\
            mean(rr_geneer_yin_cr),mean(x_FPR_turns_mean),mean(x_FDR_wufaxianlv_turns_mean),mean(y_TPR_turns_mean)]

    # df_cor_rate_record.to_csv(recom_file+"CovRecomb_"+str(seed_mut_perlin)+"_UAB"+str(len_UAB)+"_"+str(lin_diff_num)+"_result"+".csv") 
    FPR  = str(round(mean(x_FPR_turns_mean),3))
    FDR = str(round(mean(x_FDR_wufaxianlv_turns_mean),3))
    TPR = str(round(mean(y_TPR_turns_mean),3))
    return lin_diff_num,sampled_numer,mean_lin_diff,FPR,FDR,TPR
    

# Based on empirical n=1180 mutations counts from 7 countries, It is sourced from the record of a previous research (Saymon Akther, 2021).
qmatrix = {'A': { "C": 0.1083, "G": 0.7000, "T": 0.1917 },
    'C': { "A": 0.0475, "G": 0.0033, "T": 0.9492 },
    'G': { "A": 0.2102, "C": 0.0931, "T": 0.6967 },
    'T': { "A": 0.1025, "C": 0.795, "G": 0.1025 }}


seed_gen = 5 
temp_file = code_path = "/home/soniali/Desktop/03_CovRecomb/Simulation_Test/HOMO_number/"
creat_dir(temp_file)
generations_list = [6]  
turns_list = 5
sample_rate = 5
seed_gene_homo_rate = 0.1 
parallel_prop = 0.5 
seq_length = 29903
rec_rate = 0.1
mut_rate = float(8E-4/6)
len_UAB = 2
homo_num = 1

column_names = ["homo_num","seed_gen","generations_list","sampled_numer","turns_list","sample_rate","seed_gene_homo_rate","parallel_prop","seed_mut_perlin","UAB","lin_diff_num","mean_lin_diff","FPR","FDR","TPR"]
with open(temp_file+"HOMO_number.txt","a+") as f:
            f.write(",".join(column_names)+"\n")

for homo_num in [1,5,10]:
    for seed_mut_perlin in [3,5,8,10,15,20,30,40,50,60,70,80,90,100,130,150,200]:
        lin_diff_num,sampled_numer,mean_lin_diff,FPR,FDR,TPR = test_ld_genr(homo_num,seed_gen, seed_gene_homo_rate, parallel_prop,seed_mut_perlin, turns_list,generations_list,temp_file,code_path,qmatrix,sample_rate,len_UAB,seq_length,mut_rate,rec_rate)
        with open(temp_file+"HOMO_number.txt","a+") as f:
            ROC_results = [str(homo_num),str(seed_gen),str(generations_list),sampled_numer,str(turns_list),str(sample_rate),str(seed_gene_homo_rate),str(parallel_prop),str(seed_mut_perlin),str(len_UAB),str(lin_diff_num),str(mean_lin_diff),FPR,FDR,TPR]
            f.write(",".join(ROC_results)+"\n")