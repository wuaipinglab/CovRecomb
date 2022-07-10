'''
Function: Compare the elapsed time of multiple trials for the CovRecomb method and the 3SEQ method.
'''

import pandas as pd
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import os
# os.chdir(os.getcwd()+"/")


def fitting(x, y, deg, color):
    parameter = np.polyfit(x, y, deg)
    p = np.poly1d(parameter)
    plt.plot(x, p(x), color=color)


def bord_line(ax, bwith, line_color):
    ax.spines['top'].set_color(line_color)
    ax.spines['right'].set_color(line_color)
    ax.spines['left'].set_color(line_color)
    ax.spines['bottom'].set_color(line_color)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['bottom'].set_linewidth(bwith)


def obtain_sample_size(candi_file, code_path, selected_list):
    generation_sample_size = {}
    for lll in candi_file:
        sample_size = []
        for selected in selected_list:

            path = code_path + lll + "/" + "/seq_method/" + selected + "/"
            files = os.listdir(path)
            for file in files:
                if os.path.isdir(path + file):

                    for fi in os.listdir(path + file):
                        flag = 0

                        if ".3s.log" in fi:
                            flag = 1
                            parene_file_path = path + file + "/" + fi

                            with open(parene_file_path, "r") as h:
                                rows = h.readlines()
                            for row in rows:
                                if "sequences as parents." in row:
                                    pa_si = int(row.strip().split("Using ")[1].split(" sequences as parents")[0])
                                elif "sequences as children." in row:
                                    ch_si = int(row.strip().split("Using ")[1].split(" sequences as children.")[0])

                        if flag == 1:
                            sample_size.append(pa_si + ch_si)
                        else:
                            continue

        generation_sample_size[lll] = int(np.mean(sample_size))
    return generation_sample_size


def main(lin_diff_num_list, candi_file):
    bwith = 0.5
    line_color = "#c0c0c0"
    code_path = os.getcwd() + "/"
    selected_list = []
    for i in lin_diff_num_list:
        selected_list.append("lin_diff_num" + str(i))

    generation_sample_size = obtain_sample_size(candi_file, code_path, selected_list)

    # 1 Minutes = 60000000 Microseconds
    record_recomb_ld = {}
    record_seq_ld = {}
    for num in lin_diff_num_list:
        record_recom = {}
        for gen in candi_file:

            path = code_path + gen + "/covrecom/"
            gener_num = gen.split("Generation")[1]

            if "rr0.1" in path:
                gener_num = gener_num.split("_rr0.1")[0]
                cov_time = path + "ld" + str(num) + "gener" + str(gener_num) + "_truns10_result.csv"
            else:
                cov_time = path + "ld" + str(num) + "gener" + str(gener_num) + "_truns10_result.csv"

            df = pd.read_csv(cov_time)
            recom_time = df["usetime"].tolist()
            record_recom[gen] = int(mean(recom_time))
        record_recomb_ld["CovRecomb" + str(num)] = record_recom

        record_seq = {}
        temp_time = []
        for gen in candi_file:

            path = code_path + gen + "/"
            gener_num = gen.split("Generation")[1]
            seq_time = path + "seq_method/lin_diff_num" + str(num) + "/"
            files = os.listdir(seq_time)
            for file in files:
                para_record = "lin_diff_num" + str(num) + "gener" + str(gener_num)

                if os.path.exists(seq_time + file + "/" + para_record + str(file) + ".3s.log"):
                    file_path = seq_time + file + "/" + para_record + str(file) + ".3s.log"

                    with open(file_path, "r") as h:
                        rows = h.readlines()
                        for r in rows:
                            if "100.000%" in r:
                                min_second = r.strip().split("      ")[1].split("000:")[1]
                    temp_time.append(int(min_second.split(":")[0]) * 60000000 + int(min_second.split(":")[1]) * 1000000)

            record_seq[gen] = int(mean(temp_time))
        record_seq_ld["3SEQ" + str(num)] = record_seq

    Str_gene_1 = [val for val in list(generation_sample_size.keys()) for i in range(2 * len(lin_diff_num_list))]
    Str_sz_2 = [val for val in list(generation_sample_size.values()) for i in range(2 * len(lin_diff_num_list))]
    Str_method_3 = (["CovRecomb"] * len(lin_diff_num_list) + ["3SEQ"] * len(lin_diff_num_list)) * len(generation_sample_size)
    Str_ld_4 = lin_diff_num_list * len(generation_sample_size) * 2
    record_seq_ld.update(record_recomb_ld)

    Str_time_5 = []
    for i in range(len(Str_method_3)):
        elapsed_time = record_seq_ld[str(Str_method_3[i]) + str(Str_ld_4[i])][Str_gene_1[i]]
        Str_time_5.append(elapsed_time)

    df_recom_time = pd.DataFrame(columns=["generation", "sample_size", "method", "lineage_differential_mutations", "elapsed_time"])
    df_recom_time["generation"] = Str_gene_1
    df_recom_time["sample_size"] = Str_sz_2
    df_recom_time["method"] = Str_method_3
    df_recom_time["lineage_differential_mutations"] = Str_ld_4
    df_recom_time["elapsed_time"] = Str_time_5

    df_recom_time.to_csv(code_path + "/compare_time.csv")

    # plot
    your_array = array(df_recom_time["elapsed_time"])
    import numpy as np
    output = np.log(your_array)
    df_recom_time["elapsed_time_log10"] = output

    import matplotlib.pyplot as plt
    import numpy as np

    method_name = ["3SEQ", "CovRecomb"]
    for need_ld in lin_diff_num_list:
        df_recom_time_2 = df_recom_time.loc[df_recom_time["lineage_differential_mutations"] == need_ld]
        for m in method_name:
            temp_df = df_recom_time_2[df_recom_time_2["method"] == m]
            temp_df2 = temp_df[temp_df["lineage_differential_mutations"] == need_ld]
            globals()[m[1:] + "ld" + str(need_ld)] = temp_df2["elapsed_time_log10"].tolist()

        plt.cla()
        x = labels
        y1 = array(eval(str("SEQld" + str(need_ld))))
        y2 = array(eval(str("ovRecombld" + str(need_ld))))

        plt.figure(figsize=[6, 3.2], dpi=300)
        plt.scatter(x, y1, color="#DE8448")
        plt.scatter(x, y2, color="#4F75BB")

        fitting(x, y1, 2, "#DE8448")
        fitting(x, y2, 2, "#4F75BB")

        plt.title("")
        plt.xlabel('Sample size', fontweight="bold", size=10)
        plt.xticks(
            x, labels, size=10,
        )
        # max(eval(str("SEQld"+str(need_ld))))
        # min(eval(str("ovRecombld"+str(need_ld))))
        plt.ylim(ymin=10, ymax=20)
        # plt.ylim(ymin = min(min(y1),min(y1))-0.01,ymax =  max(max(y1),max(y1))+0.01)
        plt.ylabel('Elapsed time (log10 of ms)', fontweight="bold", size=10)
        plt.legend(method_name)
        ax = plt.gca()
        bord_line(ax, bwith, line_color)
        plt.savefig(code_path + "all_time" + str(need_ld) + ".pdf")


if __name__ == '__main__':
    candi_file = ["Generation4", "Generation5", "Generation8", "Generation7"]
    lin_diff_num_list = [15, 21, 38, 46]
    main(lin_diff_num_list, candi_file)
