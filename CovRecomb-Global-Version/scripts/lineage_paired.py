'''
Function: Detect the linegae or varitant paired patterns among the detected independent recombination events.

Input: Two files: Independent_recombination_event.csv, WHO_panglolin.csv.

Output: Two files: (1) lineage_paired.csv.
				   (2) variant_paired.csv.
'''

import argparse
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser(description='Read the qulified sequences.', usage='''python3 lineage_paired.py''')

    parser.add_argument("-dir", "--dirpath", default="/home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12/", help="The totoal folder address of this work.\nDefault: /home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12")
    args = parser.parse_args()
    dirpath = args.dirpath
    dirpath_pattern = dirpath + "3_recom_pattern/"

    for file in os.listdir(dirpath_pattern):
        if "Independent_recombination_event_" in file:
            df_recom = pd.read_csv(dirpath_pattern + file, header=0, delimiter=",")
            A_list = []
            B_list = []
            num = 0
            for i in df_recom.index:
                epi = df_recom.loc[i, "sample_id"]
                mode_epi = df_recom.loc[i, "mutation_pattern"]
                if mode_epi.startswith("XXXX") or mode_epi.startswith("YXXXX") or mode_epi.startswith("YYXXXX") or mode_epi.startswith("YYYXXXX"):

                    A_list.append(df_recom.loc[i, "lineage_X"])
                    B_list.append(df_recom.loc[i, "lineage_Y"])
                elif mode_epi.startswith("YYYY") or mode_epi.startswith("XYYYY") or mode_epi.startswith("XXYYYY") or mode_epi.startswith("XXXYYYY"):
                    num += 1
                    A_list.append(df_recom.loc[i, "lineage_Y"])
                    B_list.append(df_recom.loc[i, "lineage_X"])
                else:
                    print(i, epi, df_recom.loc[i, "mutation_pattern"])

    df_recom["left_lineage"] = A_list
    df_recom["right_lineage"] = B_list
    df_recom["X_Y"] = df_recom["left_lineage"] + "_" + df_recom["right_lineage"]
    ALL_MODE = df_recom["X_Y"].tolist()

    already = []
    stast_record = {}
    for c in ALL_MODE:
        if c not in already:
            count = ALL_MODE.count(c)
            already.append(c)
            stast_record[c.split("_")[0] + "_" + c.split("_")[1]] = int(count)

    # sort XY combinations according to the number of modes
    recom_num_sort = dict(sorted(stast_record.items(), key=lambda x: x[1], reverse=True))
    # for plot
    for i in recom_num_sort:
        with open(dirpath_pattern + "lineage_paired.csv", "a+") as f:
            f.write(i.split("_")[0] + "," + i.split("_")[1] + "," + str(recom_num_sort[i]) + "\n")

    # load the lineage-paired data
    df = pd.read_csv(dirpath_pattern + "lineage_paired.csv", header=None)
    df.columns = ["A", "B", "mode"]
    # change for the name of variant
    df_lin = pd.read_csv(dirpath + "0_raw_data/WHO_panglolin.csv")
    lin_suo = {}
    for i in df_lin.index:
        lin_suo[df_lin.loc[i, "pangolin"]] = df_lin.loc[i, "WHO"]

    all_lin = list(lin_suo.keys())

    mode_clus = {}
    for i in df.index:
        mode_clus["M" + str(i)] = df.loc[i, "A"], df.loc[i, "B"], df.loc[i, "mode"]

    all_AB = []
    all_AB_mode = {}
    for m in mode_clus:
        tempB = tempA = ""
        for l in all_lin:
            if l in mode_clus[m][0]:
                tempA = lin_suo[l]
            if l in mode_clus[m][1]:
                tempB = lin_suo[l]
        if tempA == "":
            tempA = mode_clus[m][0]
        if tempB == "":
            tempB = mode_clus[m][1]
        single_mode = mode_clus[m][2]
        all_AB.append((tempA, tempB))
        all_AB_mode[m] = (tempA, tempB, single_mode)
        print(tempA, tempB, single_mode)

    suolv = {}
    cout = 0
    for i in set(all_AB):
        i_m = 0
        for m in all_AB_mode:
            if (all_AB_mode[m][0], all_AB_mode[m][1]) == i:
                i_m += all_AB_mode[m][2]
        cout += i_m
        if i_m >= 2:
            suolv[i] = i_m

    df_suo = pd.DataFrame(columns=df.columns)
    df_suo["A"] = [i[0] for i in suolv.keys()]
    df_suo["B"] = [i[1] for i in suolv.keys()]
    df_suo["mode"] = [v for v in suolv.values()]
    if df_suo.shape[0] >= 1:
        df_suo.to_csv(dirpath_pattern + "variant_paired.csv", index=None)
    else:
        print(" None variant paired were identified.")


if __name__ == '__main__':
    main()
