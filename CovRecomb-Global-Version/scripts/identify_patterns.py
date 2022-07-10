'''
Function: Detect the independent recombination events from all the identified putative recombinants.

Input: Four files: snp_norm.txt, fv_norm_cluster.txt, metadata-2.tsv, cluster0.7368.csv, 

Output: Three files: (1) Independent_recombination_event.csv, the detected independent recombination events.
				     (2) Events_clustered.csv, the transmitted genomes and their corrsponding ancient-like genomes.
				     (3) mode_country_sample_id.csv, the temporal and spatial distribution and the number of recombinants for each recombination event.
'''

import pandas as pd
import copy
import datetime
import os
import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description='Read the qulified sequences.', usage='''python3 identify_patterns.py''')

    parser.add_argument("-dir", "--dirpath", default="/home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12/", help="The totoal folder address of this work.\nDefault: /home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12")
    parser.add_argument("-t", "--threshold", default=4, help="The least number of sequential feature mutations.\nDefault: 4")
    parser.add_argument("-v", "--verified_output", default="1_putative_recomb_final.csv", help="The filename of the putative recombinants already verified by the epidemiology background.\nDefault: 1_putative_recomb_final.csv\n")
    parser.add_argument("-cl", "--cluster_lineage", default=0.7368, help="The cutoff value to cluster sublineages.\nDefault: 0.7368")
    args = parser.parse_args()

    dirpath = args.dirpath
    mini_simi = args.cluster_lineage
    os.chdir(dirpath.split("data/2022_02_12/")[0] + "scripts/")
    from functions_set import creat_dir, del_star, sort_humanly, detect_mode, diff_days_fm, get_keys, extract_AB_bk, find_AB_ances, left_right_lin, get_time_epi, left_right_lin

    ############## part1  read snp for each sample
    print("\n", "----------------------- Load datasets -----------------------", "\n")
    lineage_file = dirpath + '1_filtered_data/fv_norm_cluster.txt'
    dirpath_pattern = dirpath + "3_recom_pattern/"
    creat_dir(dirpath_pattern)
    dirpath_recomb = dirpath + "2_recomb_identified/"
    verified_output = dirpath_recomb + args.verified_output
    path_breakloc = dirpath_recomb + "bk_info/"

    # read each sample's mutations
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

    print("\n", "----------------------- Find recombination events -----------------------", "\n")
    df = pd.read_csv(verified_output)
    df_recom_all = set(df['sample_id'])
    df = df.sort_values(by=["sample_id", "mutation_pattern"], ascending=("False", "False"))
    df["situation"] = df['lineage_X'] + "_" + df['lineage_Y']
    df_final = df.sort_values(by=["collection_time"], ascending=("True"))

    # build a null dataframe
    df_grouped = pd.DataFrame(columns=list(df_final))
    xxxx = 0
    each_group = {}
    # judge the parental lineages
    for i in set(df_final['situation']):
        xxxx += 1
        print("parental lineags:", i, "No.", str(xxxx))
        temp_ABC = (df_final.loc[df_final["situation"] == i]).sort_values(by=["collection_time"], ascending=("True"))
        df_first = pd.DataFrame(columns=list(df_final))
        df_candi = pd.DataFrame(columns=list(df_final))
        df_candi_in = pd.DataFrame(columns=list(df_final))

        if temp_ABC.shape[0] == 1:
            # The first sampled sequence ---- Rule 1
            df_first = pd.concat([df_first, temp_ABC])
            df_grouped = pd.concat([df_grouped, df_first])

            group_epi = list(temp_ABC["sample_id"])[0]
            group_name = i + str("_group0_") + group_epi + "_" + country[group_epi] + "_" + collected_time[group_epi]
            each_group[group_name] = [group_epi]

        else:
            df_first = pd.concat([df_first, temp_ABC.head(1)])

            #  With diverse mutations than previous genomes  ---- Rule 2
            ABC_index = (temp_ABC.index).tolist()
            left_A = del_star(temp_ABC.iloc[0]["lineage_X"])
            right_B = del_star(temp_ABC.iloc[0]["lineage_Y"])
            epi = temp_ABC.iloc[0]["sample_id"]

            first_mode = detect_mode(left_A, right_B, variants_all, epi, Lineage_v)
            turn_mode = [first_mode]
            # for the remainings
            for n in ABC_index[1:]:
                temp_mode = detect_mode(left_A, right_B, variants_all, temp_ABC.loc[n]["sample_id"], Lineage_v)
                turn_mode2 = copy.deepcopy(turn_mode)
                former_time = temp_ABC.loc[ABC_index[ABC_index.index(n) - 1]]["collection_time"]
                temp_time = temp_ABC.loc[n]["collection_time"]
                if (temp_time <= (datetime.datetime.strptime(former_time, "%Y-%m-%d") + datetime.timedelta(days=+30)).strftime("%Y-%m-%d")):
                    turn_mode, df_candi = diff_days_fm(turn_mode2, df_candi, temp_mode, temp_ABC, n, turn_mode, 4)
                else:
                    df_candi = df_candi.append(temp_ABC.loc[n])
                    turn_mode.append(temp_mode)

            # should have direct parental genomes---- Rule 3
            indi_num = 0
            epi_ances_record = {}

            for cd in df_candi.index:
                epi = df_candi.loc[cd, "sample_id"]
                left_lin, right_lin = left_right_lin(cd, df_candi)
                epiV = sort_humanly(variants_all[epi])

                for file in os.listdir(path_breakloc):
                    if epi in file:
                        small_1, big_1, small_2, big_2, bk_flag, small, big = extract_AB_bk(path_breakloc, file)

                if bk_flag == 2:
                    left_mutation = []
                    right_mutation = []
                    for v in epiV:
                        loc = int(v.split("_")[0])
                        if loc <= small_1 or loc >= big_2:
                            left_mutation.append(v)
                        elif big_1 <= loc <= small_2:
                            right_mutation.append(v)

                elif bk_flag == 1:
                    left_mutation = []
                    right_mutation = []
                    for v in epiV:
                        loc = int(v.split("_")[0])
                        if loc <= small:
                            left_mutation.append(v)
                        elif loc >= big:
                            right_mutation.append(v)

                time_epi = collected_time[epi]
                out_date = (datetime.datetime.strptime(time_epi, "%Y-%m-%d") + datetime.timedelta(days=-30)).strftime("%Y-%m-%d")
                candidate_par = (set(get_keys(country, country[epi])) & set(get_time_epi(collected_time, out_date, time_epi))) - set(df_recom_all)
                A_flag, A_ances = find_AB_ances(candidate_par, left_mutation, variants_all, pango_lineage, left_lin)
                B_flag, B_ances = find_AB_ances(candidate_par, right_mutation, variants_all, pango_lineage, right_lin)
                if A_flag == 1 and B_flag == 1:
                    df_candi_in = pd.concat([df_candi_in, df_candi.loc[df_candi.index == cd]])
                    indi_num += 1
                    epi_ances_record[epi] = ("A_ances:", A_ances, "B_ances:", B_ances)
                else:
                    continue

            # cluster samples to the identified events
            df_situation = pd.concat([df_first, df_candi_in])
            df_grouped = pd.concat([df_grouped, df_situation])

            all_epi = list(temp_ABC.sample_id)
            if df_situation.shape[0] == 1:
                group_epi = list(df_first["sample_id"])[0]
                group_name = i + str("_group0_") + group_epi + "_" + country[group_epi] + "_" + collected_time[group_epi]
                each_group[group_name] = all_epi
            else:
                recom_event = list(df_situation["sample_id"])
                recom_event_mut = {}
                for ev in recom_event:
                    recom_event_mut[ev] = sort_humanly(variants_all[ev])

                already_ev = []
                for temp_id in all_epi:
                    # judge which event should the recombinants belongs to
                    if temp_id in recom_event:
                        already_ev.append(temp_id)
                        group_name = i + "_" + temp_id + "_" + country[temp_id] + "_" + collected_time[temp_id]
                        each_group[group_name] = [temp_id]
                    else:
                        # make a loop for each event
                        epi_mut = sort_humanly(variants_all[temp_id])
                        inher_numt_num = {}
                        for cand in already_ev:
                            inher_numt_num[cand] = len(set(recom_event_mut[cand]) & set(epi_mut))

                        res = []
                        for key, value in inher_numt_num.items():
                            if value == max(inher_numt_num.values()):
                                res.append(key)

                        if len(res) == 1:
                            group_name = i + "_" + res[0] + "_" + country[res[0]] + "_" + collected_time[res[0]]
                            each_group[group_name].extend([temp_id])
                        else:
                            ca_num = {}
                            for ca in res:
                                ca_num[ca] = len(set([region[ca], country[ca], division[ca]]) & set([region[temp_id], country[temp_id], division[temp_id]]))

                            cas = []
                            for key, value in ca_num.items():
                                if value == max(ca_num.values()):
                                    cas.append(key)

                            if len(cas) == 1:
                                group_name = i + "_" + cas[0] + "_" + country[cas[0]] + "_" + collected_time[cas[0]]
                                each_group[group_name].extend([temp_id])
                            else:
                                group_name = i + "_" + cas[-1] + "_" + country[cas[-1]] + "_" + collected_time[cas[-1]]
                                each_group[group_name].extend([temp_id])

    each_group_sort = dict(sorted(each_group.items(), key=lambda item: len(item[1]), reverse=True))
    num = 0
    for n in each_group_sort:
        num += len(each_group_sort[n])

    print("Clustered recombinants:", num)
    print("The overall independent recombination event:", str(len(each_group_sort)))

    df_grouped_sort = df_grouped.sort_values(by=["collection_time"], ascending=("True"))
    df_grouped_sort.to_csv(dirpath_pattern + "Independent_recombination_event_" + str(df_grouped_sort.shape[0]) + ".csv", index=None)

    output_file_name = dirpath_pattern + "Events_clustered_" + str(len(each_group_sort)) + ".csv"
    with open(output_file_name, "a+") as f:
        f.write("Group,Recombination_event,Collection_data,Region,Country,Division,Number_of_epidemic_recombinant,Geographical_distribution,Group_members" + "\n")

    recom_event_id = list(df_grouped_sort.sample_id)
    for n in recom_event_id:
        for p in each_group_sort:
            if n not in p:
                continue
            else:
                epi_country = []
                for id in each_group_sort[p]:
                    epi_country.append(country[id])
                stast_country = {}
                for c in epi_country:
                    stast_country[c] = int(epi_country.count(c))

                geo_dis = ""
                for a in stast_country:
                    geo_dis += a + "_" + str(stast_country[a]) + "/"
                with open(output_file_name, "a+") as f:
                    f.write("_".join(p.split("_")[0:2]) + "," + n + "," + collected_time[n] + "," + region[n] + "," + country[n] + "," + division[n] + "," + str(len(each_group_sort[p])) + "," + geo_dis + "," + "/".join(each_group_sort[p]) + "\n")

    # sort the events by the number of their offspring
    each_group_sort = dict(sorted(each_group.items(), key=lambda item: len(item[1]), reverse=True))
    print("\n", "----------------------- Output recombination events -----------------------", "\n")
    with open(dirpath_pattern + "/mode_country_sample_id_" + str(len(each_group)) + ".csv", "a+") as f:
        for g in each_group_sort:
            f.write(g + "," + str(len(each_group_sort[g])) + ",")

            country_list = []
            for e in each_group_sort[g]:
                country_list.append(country[e])
            country_num = {}
            for c in set(country_list):
                country_num[c] = country_list.count(c)
            country_num_sort = dict(sorted(country_num.items(), key=lambda item: item[1], reverse=True))

            for c in country_num_sort:
                f.write(c + ":" + str(country_num_sort[c]) + ",")
            f.write("\n")


if __name__ == '__main__':
    main()