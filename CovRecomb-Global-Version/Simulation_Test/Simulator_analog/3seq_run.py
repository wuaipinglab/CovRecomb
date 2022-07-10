'''
Author: your name
Date: 2022-04-10 16:29:49
LastEditTime: 2022-07-07 21:48:28
LastEditors: Sonia-Ljy lijysunny@sina.com
Description: 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
FilePath: /undefined/home/soniali/Desktop/03_recom_0308/10_verification/lineage_zhenshi_diff/seq3_run.py
'''

import os
import subprocess
import argparse
import sys
os.chdir(os.getcwd() + "/")

def get_parent_name(row_info, loc):
    if "gen" not in row_info[loc]:
        PA = row_info[loc]
    elif "gen" in row_info[loc]:
        PA = row_info[loc].split("_")[-1]
    return PA


def main(sysargs=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        description='Simulator_CovRecombTest',
        usage='''python3 3seq_run.py -ld [15,21,38,46] -gen [4]''')

    parser.add_argument(
        "-ld",
        help= "The list number of differential mutations between simulated lineages, such as [46]",
        default=[46])
    parser.add_argument(
        "-gen",
        "--generation",
        help="The list generation number for each simulation dataset \n Default: [4]",
        default=[4])
    parser.add_argument("-cp",
                        "--code_path",
                        default=os.getcwd() + "/",
                        help="The code file")

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    lin_diff_num_list = eval(args.ld)
    generations_list = eval(args.generation)
    code_path = args.code_path
    temp_file = code_path + "compare_3seq/"
    fasta_file = temp_file + "seq_method/"

    selected_list = []
    for i in lin_diff_num_list:
        selected_list.append("lin_diff_num" + str(i))

    for g in generations_list:
        gener_num = g

    for selected in selected_list:
        num = selected.split("num")[1]
        path = fasta_file + selected + "/"
        os.chdir(path)
        files = os.listdir(path)

        print("\n", "------------------------ Run the 3SEQ method to analysis the simulated datasets ------------------------", "\n")
        for file in files:
            if os.path.isdir(path + file):
                for fi in os.listdir(path + file):
                    if "parental_seq.fasta" in fi:
                        parene_file_path = path + file + "/" + fi
                    elif "recom_seq.fasta" in fi:
                        child_file_path = path + file + "/" + fi
                        record_para = path + file + "/" + fi.split(
                            "_recom_seq")[0]

                subprocess.call([
                    "cd /home/soniali/Downloads/3seq_build/ && \
                    echo y | ./3seq -f %s %s -id %s -t1" %
                    (parene_file_path, child_file_path, record_para)
                ],
                                shell=True)

        print("\n", "------------------------ Record the coverage rate of 3SEQ method ------------------------", "\n")
        each_condi_correct_rate = {}
        for file in files:
            if os.path.isdir(path + file):
                para_record = selected + "gener" + str(
                    gener_num) + "turns" + file.split("turns")[1]

                if os.path.exists(path + file + "/" + para_record + ".3s.rec"):
                    file_path = path + file + "/" + para_record + ".3s.rec"
                elif os.path.exists(path + file + "/" + para_record +
                                    "_.3s.rec"):
                    file_path = path + file + "/" + para_record + "_.3s.rec"
                else:
                    continue

                with open(file_path, "r") as hseq:
                    rows = hseq.readlines()
                    ###
                    if len(rows) == 1:  # Not identified any recombinant
                        from Bio import SeqIO
                        Prental_ID_name = []
                        Child_ID_name = []

                        with open(
                                path + file + "/" + para_record +
                                "_parental_seq.fasta", "r") as pf:
                            for record in SeqIO.parse(pf, "fasta"):
                                Prental_ID_name.append(
                                    str(record.id).split("_")[-1])

                        with open(
                                path + file + "/" + para_record +
                                "_recom_seq.fasta", "r") as pf:
                            for record in SeqIO.parse(pf, "fasta"):
                                Child_ID_name.append(
                                    (str(record.id).split("bk")[0].split("_")
                                     [-1], str(record.id).split("bk")[1].split(
                                         "_")[-1].split("/")[0]))

                        count_recomb = 0
                        for ch in Child_ID_name:
                            if ch[0] in Prental_ID_name and ch[
                                    1] in Prental_ID_name:
                                count_recomb += 1

                        if count_recomb >= 1:
                            overall_correct_rate = 0
                        else:
                            overall_correct_rate = "Unknown"

                    else:
                        seq3_true_pos = 0
                        seq3_false_pos = 0
                        seq3_true_neg = 0
                        seq3_false_neg = 0

                        non_sig = 0
                        correct_num = 0
                        false_num = 0

                        for row in range(1, len(rows)):
                            row_info = rows[row].split("\t")
                            PA = get_parent_name(row_info, 0)
                            PB = get_parent_name(row_info, 1)

                            CH = rows[row].split("\t")[2]
                            pvalue = float(rows[row].split("\t")[9])

                            if pvalue >= 0.05:
                                non_sig += 1
                            elif pvalue < 0.05:
                                CH1 = CH.split("bk")[0].split("_")[-1]
                                if "_" in CH.split("bk")[1]:
                                    CH2 = CH.split("bk")[1].split(
                                        "/")[0].split("_")[-1]
                                else:
                                    CH2 = "lin" + CH.split("bk")[1].split(
                                        "/")[0].split("lin")[-1]

                                CH_set = {CH1, CH2}
                                PARENT_set = {PA, PB}
                                if (CH_set - PARENT_set) == set():  # CORRECT
                                    correct_num += 1
                                    if CH1 != CH2:
                                        seq3_true_pos += 1  # TP
                                    elif CH1 == CH2:
                                        seq3_true_neg += 1  # TN
                                elif (CH_set - PARENT_set) != set():  # FALSE
                                    false_num += 1
                                    if CH1 != CH2:
                                        seq3_false_pos += 1  # FP
                                    elif CH1 == CH2:
                                        seq3_false_neg += 1  # FN

                        if (false_num + correct_num + non_sig) == (len(rows) -
                                                                   1):
                            overall_correct_rate = round(
                                seq3_true_pos / (len(rows) - 1), 4)
                        else:
                            overall_correct_rate = ""

                each_condi_correct_rate[file] = overall_correct_rate

        print("\n", "------------------------ Output the coverage rate of 3SEQ in each independent trials ------------------------", "\n")
        with open(temp_file + "3seq_result_ld" + str(num) + ".txt", "w+") as seq_file:
            for k in each_condi_correct_rate:
                seq_file.write(k + ":\t" + str(each_condi_correct_rate[k]) +
                               "\n")


if __name__ == "__main__":
    main()
    