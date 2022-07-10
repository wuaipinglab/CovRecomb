'''
Function: Measure the existing probability for different sequential number of de novo mutations (sample’s mutations apart from feature mutations), and thus get a conservative threshold of the number of sequential feature mutations.

Input: Three files: the output file from STEP3: fv_norm_cluster.txt; snp_norm.txt; metadata-2.tsv

Output: A figure (denovo_chance0.7368.pdf) denotes the existing proporbility for different number of sequential feature mutations.
'''

import csv
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import matplotlib
import argparse

matplotlib.use('Agg')


def main():
    # command line interface
    parser = argparse.ArgumentParser(description='Read the qulified sequences.', usage='''python3 least_number_of_FV.py''')

    parser.add_argument("-dir", "--dirpath", default="/home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12/", help="The totoal folder address of this work.\nDefault: /home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12")
    parser.add_argument("-cl", "--cluster_lineage", default=0.7368, help="The cutoff value to cluster sublineages.\nDefault: 0.7368")

    args = parser.parse_args()

    # read feature mutations for X pool and Y pool
    dirpath = args.dirpath

    lineage_file = dirpath + "1_filtered_data/fv_norm_cluster.txt"
    mini_simi = args.cluster_lineage
    print("\n", "------------------- STRAT ------------------------------------------------", "\n")
    # read feature SNP for A pool and B pool
    Lineage_v = {}
    linA_list = []
    with open(lineage_file, 'r') as f:
        for i in csv.reader(f):
            if len(i[2:]) >= 1:
                linA_list.append(i[0])
                Lineage_v[i[0]] = i[2:]

    # read each cluster's lineages
    clu_lin = {}
    with open(dirpath + "1_filtered_data/cluster" + str(mini_simi) + ".csv", 'r') as h:
        for i in csv.reader(h):
            clu = i[0]
            lins = i[1:-1]
            clu_lin[clu] = lins

    # part2  read each sample's mutations
    variants_all = {}
    with open(dirpath + "1_filtered_data/snp_norm.txt", "r") as f:
        for i in f.readlines():
            i_2 = i.split(':')
            variants_all[i_2[0]] = i_2[1].strip().split(',')

    ########## part3  read each sample's epi_info
    pango_lineage = {}
    with open(dirpath + "0_raw_data/metadata-2.tsv", "r") as f:
        next(f)
        for i in csv.reader(f, delimiter='\t'):
            pango_lineage[i[2]] = i[18]

    # calculate the number of sequential denovo mutations in each sample
    raw_count = {'0': 0, '1': 0}
    for lin in linA_list:
        if "cluster" in lin:
            clus = clu_lin[lin.split("_")[0]]
            lin_epi = [k for k, v in pango_lineage.items() if v in clus]
        else:
            lin_epi = [k for k, v in pango_lineage.items() if v == lin]
        for i in lin_epi:
            epi_record = []
            num_denovo = 0
            try:
                for v in variants_all[i]:
                    if v in set(Lineage_v[lin]):
                        epi_record.append(num_denovo)
                        num_denovo = 0
                    else:
                        num_denovo += 1
                epi_record.append(num_denovo)
                denovo_count = dict(Counter(epi_record))
                raw_count = dict(Counter(raw_count) + Counter(denovo_count))
            except KeyError:
                continue

    lin_denovo_count = dict([(k, raw_count[k]) for k in sorted(raw_count.keys())])
    print("------------------- Chances for different threshold ------------------- " + "\n")
    print(lin_denovo_count)

    # the exsiting chance for different sequential denovo mutations
    num = 0
    for d in lin_denovo_count:
        if d != 0:
            num = num + lin_denovo_count[d]
        else:
            continue

    # the exsiting proporbility for different sequential denovo mutations
    lin_denovo_prop = {}
    for d in lin_denovo_count:
        if d != 0:
            lin_denovo_prop[d] = lin_denovo_count[d] / num
        else:
            continue

    Proportion = list(lin_denovo_prop.keys())[0:13]
    CF = list(lin_denovo_prop.values())[0:13]
    print("------------------- Proportion for different threshold ------------------- " + "\n")
    print("Number of sequential mutations:", Proportion)
    print("Existing_proporbility:", CF)

    plt.cla()
    plt.plot(Proportion, CF, color='#8696a7', marker='o')
    plt.xlabel('Number of sequential feature mutations', fontsize=10)
    plt.ylabel('Existing proporbility', fontsize=10)
    ax = plt.gca()
    x_major_locator = MultipleLocator(1)
    ax.xaxis.set_major_locator(x_major_locator)

    y_major_locator = MultipleLocator(0.05)
    ax.yaxis.set_major_locator(y_major_locator)
    plt.grid(True)
    plt.savefig(dirpath + "1_filtered_data/" + "denovo_chance" + str(mini_simi) + ".pdf", bbox_inches='tight')


if __name__ == '__main__':
    main()
