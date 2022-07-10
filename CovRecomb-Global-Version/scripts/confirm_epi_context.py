'''
Function: Confirm whether the parental lineages were pandemic one month before the candidate putative recombinants.

Input: Four files: snp_norm.txt, fv_norm_cluster.txt, metadata-2.tsv, cluster0.7368.csv, 

Output: Two files: (1) 1_putative_recomb_final.csv, the putative recombinants with epidemiology support, that is, the final results.
				   (2) 2_month_region.csv,the spatiotemporal distribution of the final identified putative recombinants.
'''

import pandas as pd
import csv
import datetime
import argparse
import os


def main():

    parser = argparse.ArgumentParser(description='Read the qulified sequences.', usage='''python3 confirm_epi_context.py''')

    parser.add_argument("-dir", "--dirpath", default="/home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12/", help="The totoal folder address of this work.\nDefault: /home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12")
    parser.add_argument("-o", "--output_file", default="0_putative_recombinants.csv", help="The filename of the output file with identified putative recombinants and their parental sequences.\nDefault: 0_putative_recombinants.csv\n")
    parser.add_argument("-v", "--verified_output", default="1_putative_recomb_final.csv", help="The filename of the putative recombinants already verified by the epidemiology background.\nDefault: 1_putative_recomb_final.csv\n")
    parser.add_argument("-cl", "--cluster_lineage", default=0.7368, help="The cutoff value to cluster lineages.\nDefault: 0.7368")
    parser.add_argument("-d", "--parental_days", default=30, help="The interval days for parental lineages.\nDefault: 30")
    args = parser.parse_args()

    dirpath = args.dirpath
    os.chdir(dirpath.split("data/2022_02_12/")[0] + "scripts/")
    from functions_set import del_star, same_cluster, get_time_epi, get_keys, creat_dir, find_bk

    dirpath_recomb = dirpath + "2_recomb_identified/"
    output_file = dirpath_recomb + args.output_file
    verified_output = dirpath_recomb + args.verified_output
    mini_simi = args.cluster_lineage
    parental_days = int(args.parental_days)
    path_breakloc = dirpath_recomb + "bk_info/"
    creat_dir(path_breakloc)
    lineage_file = dirpath + '1_filtered_data/fv_norm_cluster.txt'

    # read each sample's mutations
    print("\n", "----------------------- Load datasets -----------------------", "\n")
    Strain_list_snp = []
    variants_all = {}
    with open(dirpath + "1_filtered_data/snp_norm.txt", "r") as f:
        for i in f.readlines():
            i_2 = i.split(':')
            Strain_list_snp.append(i_2[0])
            variants_all[i_2[0]] = i_2[1].strip().split(',')

    # read each sample's epi_info
    Strain_list_meta = []
    collected_time = {}
    region = {}
    country = {}
    division = {}
    pango_lineage = {}
    with open(dirpath + "0_raw_data/metadata-2.tsv", "r") as f:
        next(f)
        for i in csv.reader(f, delimiter='\t'):
            Strain_list_meta.append(i[2])
            collected_time[i[2]] = i[4]
            region[i[2]] = i[9]
            country[i[2]] = i[10]
            division[i[2]] = i[11]
            pango_lineage[i[2]] = i[18]

    cluster_dict = {}
    with open(dirpath + "1_filtered_data/cluster" + str(mini_simi) + ".csv", "r") as f:
        cluster = csv.reader(f)
        for i in cluster:
            cluster_dict[i[0]] = i[1:-1]

    Lineage_v = {}
    with open(lineage_file, 'r') as f:
        for i in csv.reader(f):
            if "cluster" not in i[0]:
                Lineage_v[i[0]] = i[2:]
            else:
                cluster_n = i[0].split("_")[0]
                for c in cluster_dict[cluster_n]:
                    Lineage_v[c] = i[2:]

    cluster_lin = {}
    with open(dirpath + "1_filtered_data/cluster" + str(mini_simi) + ".csv", "r") as f:
        for i in f.readlines():
            i_2 = i.strip().split(',')
            cluster_lin[i_2[0]] = set(i_2[1:-1])

    # 判断lineage各自属于哪个cluster中
    lin_cluster = {}
    for c in cluster_lin:
        for l in cluster_lin[c]:
            lin_cluster[l] = c

    print("\n", "----------------------- Check epidemiology background -----------------------", "\n")
    df = pd.read_csv(output_file, header=0)
    df_sortID = df.sort_values(by=['sample_id', 'mutation_pattern', "lineage_X", "lineage_Y"], axis=0, ascending=[True, True, True, True], inplace=False)

    df_recom_all = df_sortID["sample_id"].tolist()
    df_recom_A = df_sortID["lineage_X"].tolist()
    df_recom_B = df_sortID["lineage_Y"].tolist()

    count = 0
    count_epiinfo_in = 0
    col_names = ['sample_id', 'collection_time', 'region', 'country', 'division', 'pangolin', 'lineage_X', 'lineage_Y']
    df_epi_info_in = pd.DataFrame(columns=col_names)
    for i in range(df.shape[0]):
        count += 1
        print(count)
        epi_name = df_recom_all[i]

        time_epi = collected_time[epi_name]
        country_epi = country[epi_name]

        epi_A = []
        epi_B = []
        epi_A.append(del_star(df_recom_A[i]))
        epi_B.append(del_star(df_recom_B[i]))
        epi_A = same_cluster(epi_A, lin_cluster, cluster_lin)
        epi_B = same_cluster(epi_B, lin_cluster, cluster_lin)

        A_num = 0
        B_num = 0
        out_date = (datetime.datetime.strptime(time_epi, "%Y-%m-%d") + datetime.timedelta(days=-parental_days)).strftime("%Y-%m-%d")

        candidate_par = (set(get_keys(country, country_epi)) & set(get_time_epi(collected_time, out_date, time_epi))) - set(df_recom_all)
        for n in candidate_par:
            if A_num == 0 or B_num == 0:
                if pango_lineage[n] in epi_A:
                    A_num += 1
                if pango_lineage[n] in epi_B:
                    B_num += 1
            else:
                break

        if A_num == 0 or B_num == 0:
            count_epiinfo_in += 1
        else:
            df_epi_info_in = df_epi_info_in.append(df_sortID[df_sortID["sample_id"] == epi_name])
            df_epi_info_in["collection_time"] = collected_time[epi_name]
            df_epi_info_in["region"] = region[epi_name]
            df_epi_info_in["country"] = country[epi_name]
            df_epi_info_in["division"] = division[epi_name]
            df_epi_info_in["pangolin"] = pango_lineage[epi_name]
            continue

    df_epi_info_in.to_csv(verified_output, sep=',', index=False)

    print("\n", "------------------ Output the breakpoint of each putative recombinant ---------------------", "\n")
    for epi in list(df_epi_info_in["sample_id"]):
        find_bk(df_epi_info_in, epi, variants_all, path_breakloc)

    print("\n", "------------------ Analysis the spatio-temporal distribution ---------------------", "\n")
    UA_time = df_epi_info_in.sort_values(by="collection_time", ascending=True)
    region_list = UA_time["region"]
    month_list = []
    time_list = UA_time["collection_time"].tolist()
    for i in time_list:
        j = i[0:7]
        month_list.append(j)

    UA_time["time_month"] = month_list
    cut_time = ["2020-01","2020-02","2020-03","2020-04","2020-05","2020-06","2020-07","2020-08","2020-09","2020-10","2020-11","2020-12",\
        "2021-01","2021-02","2021-03","2021-04","2021-05","2021-06","2021-07","2021-08","2021-09","2021-10","2021-11","2021-12",\
        "2022-01","2022-02"]

    need_time = cut_time[cut_time.index(month_list[0]):cut_time.index(month_list[-1]) + 1]
    nt_r_all = {}
    for nt in need_time:
        UA_time_need = UA_time[UA_time["time_month"] == nt]
        nt_region = set(UA_time_need["region"].tolist())

        r_all = {}
        for nt_r in nt_region:
            nt_r_num = 0
            for j in UA_time_need["region"]:
                if j == nt_r:
                    nt_r_num += 1
            r_all[nt_r] = nt_r_num
        nt_r_all[nt] = r_all

    df_barplot = pd.DataFrame(columns=list(set(region_list)), index=need_time)
    for t in nt_r_all:
        df_barplot.loc[t] = pd.Series(nt_r_all[t])
    df_barplot.fillna(0, inplace=True)
    df_barplot = df_barplot.astype(int)

    sum_re = 0
    for i in df_barplot:
        for j in df_barplot.index:
            sum_re += int(df_barplot.loc[j, i])

    df_barplot.to_csv(dirpath_recomb + "2_month_region_" + str(df_epi_info_in.shape[0]) + ".csv", sep=',', index=True)


if __name__ == '__main__':
    main()
