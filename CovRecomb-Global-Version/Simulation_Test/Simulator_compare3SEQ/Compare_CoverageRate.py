'''
Author: Jiaying Li
Date: 2022-04-10 19:40:10
LastEditTime: 2022-10-06 16:30:25
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


def get_index(lst=None, item=''):
     return [i for i in range(len(lst)) if lst[i] == item]


def del_indel(seq3,index_list,need_bar,raw_len):
    select_n = []
    del_count = 0
    for n in range(len(seq3)):
        if n in index_list and del_count < int(raw_len-need_bar):
            del_count+=1
            continue
        else:
            select_n.append(seq3[n])
            
    return select_n


def set_label(rects):
    for rect in rects:
        height = float(format(rect.get_height(),".2f"))
        plt.text(x = rect.get_x() + rect.get_width()/2,
                y = height + 0.05, 
                s = height, 
                ha = 'center',size = 5) 


def bord_line(ax,bwith, line_color):
    ax.spines['top'].set_color(line_color)
    ax.spines['right'].set_color(line_color)
    ax.spines['left'].set_color(line_color)
    ax.spines['bottom'].set_color(line_color)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['bottom'].set_linewidth(bwith)
    
    
def draw_compare_seq_recom(num,ld_diff, temp_file, seq_file, recom_file,need_bar,line_color,bwith):
    width = 0.2
    x0_all = ["rep1","rep2","rep3","rep4","rep5","rep6"]
    x0 = x0_all[0:need_bar]
    name_list = []
    for i in x0:
        name_list.append(i.split("rep")[1])
    seq3_0 = {}
    with open(temp_file+seq_file,"r") as f:
        rows = f.readlines()
        for i in rows:
            seq3_0[i.strip().split(":")[0]] =float(i.strip().split(":")[1].split("\t")[1])

    dd = pd.read_csv(temp_file+recom_file)
    cov_recom_0 = {}
    for i in range(len(dd["Unnamed: 0"].tolist())):
        cov_recom_0[dd["Unnamed: 0"].tolist()[i].split("gene"+str(num))[1].split("RT")[0]] =  dd["y_TPR"].tolist()[i]


    ## Get the "need_bar" higest value of 3SEQ and CovRecomb
    sum_seq3_Recomb = {}
    value_list = []
    for i in cov_recom_0:
        sum_seq3_Recomb[i] = float(cov_recom_0[i])+float(seq3_0[i])
        value_list.append(float(cov_recom_0[i])+float(seq3_0[i]))

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
    
    plt.cla()
    plt.figure(figsize = [4.2,3.2],dpi = 300)
    import numpy as np
    x = np.arange(len(seq3))
    rects1 = plt.bar(x - width/2, seq3, width,color = "#DE8448") 
    rects2 = plt.bar(x + width/2, cov_recomb, width,color = "#4F75BB")
    
    # Title, legend
    plt.xticks(x,name_list)
    font = {'weight': 'normal', 'size': 6}
    plt.legend(['3SEQ','CovRecomb'],prop=font)
    # Annations
    plt.xlabel('',size = 10,fontweight = "bold")
    plt.ylim(ymin = 0,ymax = 1.12)
    plt.ylabel('Coverage rate',size = 9,fontweight = "bold")
   
    y = MultipleLocator(0.2)  
    ax = plt.gca()
    ax.yaxis.set_major_locator(y)
    bord_line(ax,bwith, line_color)
    plt.savefig(temp_file+"compare_"+str(ld_diff)+".pdf")


def main(sysargs = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='Simulator_CovRecombTest',
    # usage='''nohup CovRecomb mode_222_mafft.fasta > CovRecomb_record.log 2>&1 &''')
    usage = '''python3 Compare_CoverageRate.py -gen [4,6]''')

    parser.add_argument("-gen", help="The generation number for each simulation dataset in the format of a list \n Default: [4,6]", default=[4,6])
    parser.add_argument("-top","--need_bar", default = 6, help="The top value from turns number of trials")
    parser.add_argument("-cp","--code_path", default = os.getcwd()+"/", help="The code file")
    
    if len(sysargs)<1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    
    code_path = args.code_path
    temp_file = code_path+"CovRecomb_compare/"
    need_bar = int(args.need_bar)
    generations_list = eval(args.gen)
    
    bwith = 0.5 
    line_color = "#c0c0c0"
    for num in generations_list:
        ld_diff = "gen"+str(num)
        seq_file = "3seq_result_gener"+str(num)+".txt"
        recom_file = "seed_mut_perlin_8_gener"+str(num)+"_truns10_result.csv"
        draw_compare_seq_recom(num,ld_diff, temp_file, seq_file, recom_file,need_bar,line_color,bwith)


if __name__ == "__main__":
    main()
    