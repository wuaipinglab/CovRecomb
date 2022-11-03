'''
The core algorithm for CovRecomb method
'''

from function_set import mutation_lin_unique


def unique_lin(aftertime_lin, Lineage_v, U_mutate_unique):
    '''
    Check whether the sequential mutations were unique in the candidate parental lineages (or its descendants or siblings).

    :param aftertime_lin: the candidate parental lineage
    :param Lineage_v: dict, key: the candidate parental lineage; value: the feature mutations
    :param U_mutate_unique: list, the feature mutations combination in one of the candidate parental lineage
    
    :return: int, the flag representing the mutations combination is unique or not 
    '''
    num = 0
    same_ancient = []
    for can_A in aftertime_lin:
        if set(Lineage_v[can_A]) >= set(U_mutate_unique):
            same_ancient.append(can_A)
            num += 1

    if num == 1:
        return num
    elif num > 1:
        same_A_num = mutation_lin_unique(same_ancient)
        if same_A_num == len(same_ancient):
            num = 1
        else:
            num = 2

    return num


def calcul_bk(lin_A_draw, lin_B_draw, Lineage_v, epiV):
    '''
    Check whether the sequential mutations were unique in the candidate parental lineages (or its descendants or siblings).

    :param lin_A_draw: the candidate parental lineage X
    :param lin_B_draw: the candidate parental lineage Y
    :param Lineage_v: dict, key: the candidate parental lineage; value: the feature mutations
    :param epiV: list, the mutations of the analyzing sequence
    
    :return: str, the sequential feature mutations pattern
    :return: list, the unique feature mutations of candidate parental lineage X
    :return: list, the unique feature mutations of candidate parental lineage Y
    :return: list, the feature mutations shared by the both candidate parental lineages (X and Y)
    :return: list, the unique mutations of the analyzing sequence 
    '''
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


def qu_cluster(lin_A):
    if "cluster" in lin_A:
        renew_name = lin_A.split("_")[1]+"*"
    else:
        renew_name = lin_A
    return renew_name


def recombination_detection(len_UAB, max_bk_num, must_inA, must_inB, linA_list_deep, Strain_list_def,\
    variants_all, feature_mutations, Lineage_v, mutaions_num, output_file):
    '''
    To detect whether the input sequences are the putative inter-lineage recombinant using the hypergeometric-test-based algorithm
    '''
    from scipy.stats import hypergeom
    from function_set import sort_humanly, min_pairs, sort_humanly, bk_count, obtain_pattern, parental_lineage

    for epi in Strain_list_def:
        epi_record = {}
        epiV = variants_all[epi]
        epi_feat = len(set(epiV) & set(feature_mutations))

        if epi_feat < 2*len_UAB:
            continue
        else:
            # P-value for Non-recombination
            pmin = 1
            for lin_A in linA_list_deep:
                all_AA = len(Lineage_v[lin_A]) * 2
                all_AA_epi = len(set(Lineage_v[lin_A]) & set(epiV))

                pVal = hypergeom.sf(all_AA_epi - 1, mutaions_num, all_AA,
                                    epi_feat)
                if pVal >= pmin:
                    continue
                else:
                    pmin = pVal
                    epi_record[str(lin_A) + "_" + str(lin_A)] = pVal

            # the least p-value for the Non-recombinant
            most_one = min_pairs(epi_record)
            pmin = most_one[0][1]

            epi_record = {}
            for mo in most_one:
                epi_record[mo[0]] = mo[1]

            # P-value for Recombinant (A+B/A+B+A)
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

                        if len(B_epi) < len_UAB or len(B_epi -
                                                       A_epi) < len_UAB or len(
                                                           A_epi -
                                                           B_epi) < len_UAB:
                            continue
                        else:
                            all_AB = len(Lineage_v[A]) + len(Lineage_v[B])
                            all_AB_epi = len((set(Lineage_v[A])
                                              | set(Lineage_v[B])) & set(epiV))

                            pVal = hypergeom.sf(all_AB_epi - 1, mutaions_num,
                                                all_AB, epi_feat)
                            if pVal > pmin:
                                continue
                            else:
                                pmin = pVal
                                unique_A = A_epi - B_epi
                                unique_B = B_epi - A_epi

                                if len(unique_A) < len_UAB or len(
                                        unique_B) < len_UAB:
                                    continue
                                else:
                                    union_AB_set = set(Lineage_v[A]) ^ set(
                                        Lineage_v[B])
                                    AB_epi = sort_humanly(
                                        list(union_AB_set & set(epiV)))
                                    recom_pattern = obtain_pattern(
                                        AB_epi, unique_A, unique_B)
                                    if (must_inA not in recom_pattern) or (
                                            must_inB not in recom_pattern):
                                        continue
                                    else:
                                        change = bk_count(recom_pattern)
                                        if change > max_bk_num:
                                            continue
                                        else:
                                            epi_record[str(A) + "_" +
                                                       str(B)] = pVal
            most_two = min_pairs(epi_record)

            if most_one == most_two:  # Mostly, this sample is non-recombinant
                continue
            else:
                epiV = sort_humanly(epiV)
                for M in most_two:
                    lin_A_draw, lin_B_draw = parental_lineage(M)
                    lin_record, UA_mutate_unique, UB_mutate_unique, shared_mut, denovo_mut = calcul_bk(
                        lin_A_draw, lin_B_draw, Lineage_v, epiV)

                    lin_A_draw, lin_B_draw = qu_cluster(lin_A_draw), qu_cluster(lin_B_draw)
                    # numA = unique_lin(linA_list_deep, Lineage_v,
                    #                   UA_mutate_unique)
                    # numB = unique_lin(linA_list_deep, Lineage_v,
                    #                   UB_mutate_unique)

                    # if numA == numB == 1:
                    with open(output_file, "a+") as file_epi:
                        file_epi.write(epi + "," + lin_A_draw + "," +
                                        lin_B_draw + "," + lin_record +
                                        "," + "/".join(UA_mutate_unique) +
                                        "," + "/".join(UB_mutate_unique) +
                                        "," + "/".join(shared_mut) + "," +
                                        "/".join(denovo_mut) + "\n")
