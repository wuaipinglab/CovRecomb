'''
Function: Cluster lineages with similar lineage-defining mutations with a cutoff value derived from the minimum mutual similarity among viral lineages in Alpha variant.

Input: One file: the output file from STEP2: fv_norm.txt.

Output: Three files: fv_norm_cluster.txt, cluster0.7368.csv, simi_group0.7368.csv. The decimal in the latter two filenames represents the least threshold (minimum mutual similarity) to cluster lineages.
'''

import csv
import pandas as pd
import numpy as np
import os
import argparse


def main():
    # command line interface
    parser = argparse.ArgumentParser(description='Read the qulified sequences.', usage='''python3 cluster_lineage.py''')

    parser.add_argument("-dir", "--dirpath", default="/home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12/", help="The totoal folder address of this work.\nDefault: /home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12")
    args = parser.parse_args()

    dirpath = args.dirpath
    fv_path = dirpath + "1_filtered_data/fv_norm.txt"

    os.chdir(dirpath.split("data/2022_02_12/")[0] + "scripts/")
    from functions_set import sort_humanly

    print("\n", "------------------- STRAT ------------------------------------------------", "\n")

    # read feature mutations
    Lineage_v = {}
    lineage_n = {}
    linA_list = []
    with open(fv_path, 'r') as f:
        for i in csv.reader(f):
            linA_list.append(i[0])
            lineage_n[i[0]] = i[1]
            Lineage_v[i[0]] = i[2:]

    # calculate the similarity of lineage X and Y
    sample_num = pd.DataFrame(index=linA_list, columns=linA_list)

    for A in linA_list:
        for B in linA_list:
            if np.isnan(sample_num.loc[B, A]) == True:
                mutat_A = set(Lineage_v[A])
                mutat_B = set(Lineage_v[B])
                at_least_similarity = len(mutat_A & mutat_B) / max(len(mutat_A), len(mutat_B))
                sample_num.loc[B, A] = sample_num.loc[A, B] = at_least_similarity

    # calculate the cutoff value of the Alpha variant
    alpha_list = ["B.1.1.7", "Q.1", "Q.2", "Q.3", "Q.4", "Q.5", "Q.6", "Q.7", "Q.8"]
    simi_list = []
    for a in alpha_list:
        for b in alpha_list:
            if a != b:
                simi_list.append(sample_num.loc[a, b])

    print("\n", "The threathold of cluster is:", min(simi_list), "\n")

    couple = []
    all_related = []
    for i in sample_num:
        for j in sample_num:
            if j != i and sample_num.loc[i, j] >= min(simi_list):
                couple.append([i, j])
                if i not in all_related:
                    all_related.append(i)
                if j not in all_related:
                    all_related.append(j)
            else:
                continue

    # for each (has line) point, extract the related points with it
    group = {}
    for l in all_related:
        each_group = [l]
        for i in couple:
            if i[1] == l:
                each_group.append(i[0])
            elif i[0] == l:
                each_group.append(i[1])
        group[l] = set(each_group)

    with open(dirpath + "1_filtered_data/simi_group" + str(round(min(simi_list), 4)) + ".csv", "a+") as f:
        for i in group:
            f.write(i + ":")
            for v in group[i]:
                f.write(v + ",")
            f.write("\n")

    already = []
    concluded_g = {}
    for g in group:
        if g not in already:
            concluded = []
            for m in group:
                if m != g and group[m] == group[g]:
                    concluded.append(m)
                    already.append(m)
                else:
                    continue
            concluded_g[g] = concluded

    count = 0
    cluster = {}
    for i in concluded_g:
        if len(concluded_g[i]) >= 1:
            count += 1
            groupname = "cluster" + str(int(count))
            raw = list(set(concluded_g[i]))
            raw.append(i)
            raw.sort()
            cluster[groupname] = raw
        else:
            continue

    print("\n", "The total number of cluster is:", count, "\n")

    # write the cluster info
    with open(dirpath + "1_filtered_data/cluster" + str(round(min(simi_list), 4)) + ".csv", "a+") as f:
        for i in cluster:
            f.write(i + ",")
            for v in cluster[i]:
                f.write(v + ",")
            f.write("\n")

    # record each lineage belongs to which cluster
    lin_CLA = {}
    for i in cluster:
        for A in cluster[i]:
            lin_CLA[A] = i

    # select the represent lineage(who with the most genomes)
    repres_lin = {}
    for clu in cluster:
        lingeage_n = {}
        with open(fv_path, 'r') as f:
            for i in csv.reader(f):
                if len(i[2:]) >= 1:
                    if i[0] in cluster[clu]:
                        lingeage_n[i[0]] = int(i[1])
        repres_lin[clu] = max(lingeage_n, key=lingeage_n.get)

    # write down the new FV_norm_cluster.txt file
    count_raw = 0
    count_cluster = 0
    cluster_already = []
    with open(fv_path, "r") as f:
        for i in csv.reader(f):
            if i[0] not in lin_CLA:
                count_raw += 1
                with open(dirpath + "1_filtered_data/fv_norm_cluster.txt", "a+") as h:
                    h.write(i[0] + "," + i[1] + "," + ",".join(i[2:]) + "\n")
            else:
                cluster_name = lin_CLA[i[0]]
                count_cluster += 1
                if cluster_name not in cluster_already:
                    cluster_already.append(cluster_name)
                    culster_rep_lin = repres_lin[cluster_name]
                    culster_num = lineage_n[culster_rep_lin]
                    culster_mutations = sort_humanly(set(Lineage_v[culster_rep_lin]))

                    with open(dirpath + "1_filtered_data/fv_norm_cluster.txt", "a+") as h:
                        h.write(cluster_name + "_" + culster_rep_lin + "," + culster_num + "," + ",".join(culster_mutations) + "\n")

    print("\n", "------------------- Cluster_lineages_completed ---------------------------", "\n")


if __name__ == '__main__':
    main()
