'''
Input: The coverage rate of identifying the simulated recombinants for CovRecomb and the 3SEQ method

Parameters: (1) "-ld": The number of differential mutations between simulated lineages in the format of a list; Default: [46].
			(2) "-gen","--generation": The generation number for each simulation dataset in the format of a list; Default: [4].
			(3) "-top","--need_bar": The top ranked value from multiple turns of trials; Default = 6.
			(9) "-cp","--code_path": The file address of the folder with scripts; Default = os.getcwd(+"/").

Function: Compare the coverage rate of multiple trials for the CovRecomb method and the 3SEQ method. 

Output: Figure(s) for the comparison between the coverage rate of CovRecomb and 3SEQ method in multiple trails.
'''

import pandas as pd
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import os
import argparse
import sys

# os.chdir(os.getcwd()+"/")


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
        plt.text(
            x=rect.get_x() + rect.get_width() / 2,
            y=height + 0.05,  # 竖直坐标
            s=height,  # ⽂本
            ha='center',
            size=5)  # ⽔平居中


def bord_line(ax, bwith, line_color):
    ax.spines['top'].set_color(line_color)
    ax.spines['right'].set_color(line_color)
    ax.spines['left'].set_color(line_color)
    ax.spines['bottom'].set_color(line_color)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['bottom'].set_linewidth(bwith)


def draw_compare_seq_recom(ld_diff, path, seq_file, recom_file, need_bar, gener_num, line_color, bwith):
    width = 0.2
    x0_all = ["rep1", "rep2", "rep3", "rep4", "rep5", "rep6"]
    x0 = x0_all[0:need_bar]
    name_list = []
    for i in x0:
        name_list.append(i.split("rep")[1])
    seq3_0 = []
    with open(path + seq_file, "r") as f:
        rows = f.readlines()
        for i in rows:
            seq3_0.append(float(i.strip().split(":")[1].split("\t")[1]))

    dd = pd.read_csv(path + "covrecom/" + recom_file)
    cov_recom_0 = dd["all_recom_correct_rate"].tolist()
    raw_len = len(seq3_0)

    ## Get the "need_bar" higest value of 3SEQ and CovRecomb
    sum_seq3_Recomb = []
    for i in range(len(cov_recom_0)):
        sum_seq3_Recomb.append(float(cov_recom_0[i]) + float(seq3_0[i]))

    import heapq
    max_number = heapq.nsmallest(raw_len - need_bar, sum_seq3_Recomb)

    max_index = []
    for t in max_number:
        index = get_index(sum_seq3_Recomb, t)
        max_index.extend(index)

    seq3 = del_indel(seq3_0, list(set(max_index)), need_bar, raw_len)
    cov_recomb = del_indel(cov_recom_0, list(set(max_index)), need_bar, raw_len)

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
    plt.xlabel('', size=10, fontweight="bold")
    plt.ylim(ymin=0, ymax=1.12)
    plt.ylabel('Coverage rate', size=9, fontweight="bold")

    y = MultipleLocator(0.2)
    ax = plt.gca()
    ax.yaxis.set_major_locator(y)
    bord_line(ax, bwith, line_color)
    plt.savefig(path + "compare_" + str(ld_diff) + "_gen" + str(gener_num) + ".pdf")


def main(sysargs=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        description='Simulator_CovRecombTest',
        usage='''python3 Compare_CoverageRate.py -ld [15,21,38,46] -gen [4] -top 6''')

    parser.add_argument("-ld", help="The number of differential mutations between simulated lineages in the format of a list \n Default: [46]", default=[46])
    parser.add_argument("-gen", "--generation", help="The generation number for each simulation dataset in the format of a list \n Default: [4]", default=[4])
    parser.add_argument("-top", "--need_bar", help="The top ranked value from multiple turns of trials", default=6)
    parser.add_argument("-cp", "--code_path", default=os.getcwd() + "/", help="The code file")

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    lin_diff_num_list = eval(args.ld)
    generations_list = eval(args.generation)
    code_path = args.code_path
    temp_file = code_path + "compare_3seq/"
    need_bar = int(args.need_bar)

    gener_num = generations_list[0]
    bwith = 0.5
    line_color = "#c0c0c0"
    for num in lin_diff_num_list:
        ld_diff = "ld" + str(num)
        seq_file = "3seq_result_ld" + str(num) + ".txt"
        recom_file = "ld" + str(num) + "gener" + str(gener_num) + "_truns10_result.csv"
        draw_compare_seq_recom(ld_diff, temp_file, seq_file, recom_file, need_bar, gener_num, line_color, bwith)

    # Rename the folder according to the parameter number of generation.
    new_filename = temp_file.split("compare_3seq/")[0] + "Generation" + str(gener_num) + "/"
    os.rename(temp_file, new_filename)


if __name__ == "__main__":
    main()
