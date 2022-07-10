'''
Author: Jiaying Li
Date: 2022-04-10 19:40:10
LastEditTime: 2022-07-07 22:08:00
LastEditors: Sonia-Ljy lijysunny@sina.com
Description: 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
FilePath: /undefined/home/soniali/Desktop/03_recom_0308/10_verification/lineage_zhenshi_diff/compare_3seq/compare_3seq_diff.py
'''

import pandas as pd
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import os
import argparse
import sys

# os.chdir(os.getcwd()+"/"
def get_index(lst=None, item=''):
    return [i for i in range(len(lst)) if lst[i] == item]


def del_indel(seq3, index_list, need_bar, raw_len):

    select_n = []
    del_count = 0
    for n in range(len(seq3)):
        if n in index_list and del_count < int(raw_len - need_bar):
            del_count += 1
            continue
        else:
            select_n.append(seq3[n])

    return select_n


def set_label(rects):
    for rect in rects:
        height = float(format(rect.get_height(), ".2f"))
        plt.text(x=rect.get_x() + rect.get_width() / 2, y=height + 0.05, s=height, ha='center', size=5)


def bord_line(ax, bwith, line_color):
    ax.spines['top'].set_color(line_color)
    ax.spines['right'].set_color(line_color)
    ax.spines['left'].set_color(line_color)
    ax.spines['bottom'].set_color(line_color)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['bottom'].set_linewidth(bwith)


def draw_compare_seq_recom(ld_diff, temp_file, seq_file, recom_file, need_bar, gener_num, line_color, bwith):
    width = 0.2
    x0_all = ["rep1", "rep2", "rep3", "rep4", "rep5", "rep6"]
    x0 = x0_all[0:need_bar]
    name_list = []
    for i in x0:
        name_list.append(i.split("rep")[1])
    seq3_0 = {}
    with open(temp_file + seq_file, "r") as f:
        rows = f.readlines()
        for i in rows:
            seq3_0[i.strip().split(":")[0]] = float(i.strip().split(":")[1].split("\t")[1])

    dd = pd.read_csv(temp_file + "covrecom/" + recom_file)
    cov_recom_0 = {}
    for i in range(len(dd["Unnamed: 0"].tolist())):
        cov_recom_0["turns" + dd["Unnamed: 0"].tolist()[i].split("turns")[1]] = dd["all_recom_correct_rate"].tolist()[i]

    sum_seq3_Recomb = {}
    value_list = []
    for i in cov_recom_0:
        sum_seq3_Recomb[i] = float(cov_recom_0[i]) + float(seq3_0[i])
        value_list.append(float(cov_recom_0[i]) + float(seq3_0[i]))

    value_list_sort = sorted(value_list)[-int(need_bar):]

    seq3 = []
    cov_recomb = []
    for i in sum_seq3_Recomb:
        if sum_seq3_Recomb[i] in value_list_sort:
            seq3.append(seq3_0[i])
            cov_recomb.append(cov_recom_0[i])

    if len(seq3) > need_bar:
        for i in sum_seq3_Recomb:
            if sum_seq3_Recomb[i] == sorted(value_list)[-int(need_bar):][0]:
                seq3.remove(seq3_0[i])
                cov_recomb.remove(cov_recom_0[i])
                break

    # plot
    plt.cla()
    plt.figure(figsize=[4.2, 3.2], dpi=300)
    import numpy as np
    x = np.arange(len(seq3))
    plt.bar(x - width / 2, seq3, width, color="#DE8448")
    plt.bar(x + width / 2, cov_recomb, width, color="#4F75BB")

    plt.xticks(x, name_list)
    font = {'weight': 'normal', 'size': 6}
    plt.legend(['3SEQ', 'CovRecomb'], prop=font)
    # Annations
    plt.xlabel('', size=10, fontweight="bold")
    plt.ylim(ymin=0, ymax=1.12)
    plt.ylabel('Coverage rate', size=9, fontweight="bold")

    y = MultipleLocator(0.2)
    ax = plt.gca()
    ax.yaxis.set_major_locator(y)
    bord_line(ax, bwith, line_color)
    plt.savefig(temp_file + "compare_" + str(ld_diff) + "_gen" + str(gener_num) + ".pdf")


def main(sysargs=sys.argv[1:]):
    parser = argparse.ArgumentParser(description='Simulator_CovRecombTest', usage='''python3 Compare_CoverageRate.py -ld [15,21,38,46] -gen [4] -top 6''')
    parser.add_argument("-ld", help="The list number of differential mutations between simulated lineages, such as [46]", default=[46])
    parser.add_argument("-gen", "--generation", help="The list generation number for each simulation dataset \n Default: [4]", default=[4])
    parser.add_argument("-top", "--need_bar", default=6, help="The top value from turns number of trials")
    parser.add_argument("-cp", "--code_path", default=os.getcwd() + "/", help="The code file")
    parser.add_argument("-sr", "--sample_rate", help="How many number of recombinants will be sampled for one sequence ", type=int, default=100)

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    lin_diff_num_list = eval(args.ld)
    generations_list = eval(args.generation)
    code_path = args.code_path
    temp_file = code_path + "compare_3seq/"
    sample_rate = int(args.sample_rate)
    need_bar = int(args.need_bar)

    gener_num = generations_list[0]
    bwith = 0.5
    line_color = "#c0c0c0"
    for num in lin_diff_num_list:
        ld_diff = "ld" + str(num)
        seq_file = "3seq_result_ld" + str(num) + ".txt"
        recom_file = "ld" + str(num) + "gener" + str(gener_num) + "_truns10_result.csv"
        # os.chdir(os.getcwd()+"/")
        draw_compare_seq_recom(ld_diff, temp_file, seq_file, recom_file, need_bar, gener_num, line_color, bwith)

    # rename the folder
    new_filename = temp_file.split("compare_3seq/")[0] + "Sampler_Generation" + str(gener_num) + "_SR_" + str(sample_rate) + "/"
    os.rename(temp_file, new_filename)


if __name__ == "__main__":
    main()
