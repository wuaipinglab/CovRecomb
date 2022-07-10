'''
Function: Maske the start (1 - 265 nt) and the end (29,674 - 29,903 nt) untranslated region (UTR) regions in idenfied mutaions.

Input: Two files: the output files from STEP1: snp.txt, fv.txt.

Output: Two files: normalized mutations for qualified samples (snp_norm.txt) and the normalized lineage-defining feature mutation (fv_norm.txt).
'''
import csv
import copy
import argparse
import os


def normalize_position(j, geno_min, geno_max):
    V_copy = copy.deepcopy(j)
    if j != "" and j != ['']:
        for v in j:
            pos = int(v.split("_")[0])
            if pos < geno_min or pos > geno_max:
                V_copy.remove(v)
    else:
        V_copy = V_copy

    return str(V_copy)


def main():
    # command line interface
    parser = argparse.ArgumentParser(description='Read the qulified sequences.', usage='''python3 pos_norm.py''')

    parser.add_argument("-dir", "--dirpath", default="/home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12/", help="The totoal folder address of this work.\nDefault: /home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12")
    parser.add_argument("-del", "--delete_lineage", default=["BA.3"], help="The lineage(s) should be deleted which might came from recombination.\nDefault: ['BA.3']")
    parser.add_argument("-t", "--threshold", default=4, help="The least number of sequential feature mutations.\nDefault: 4")
    parser.add_argument("-mmin", "--masked_min", default=266, help="The minimum position for site need to be masked.\nDefault: 266")
    parser.add_argument("-mmax", "--masked_max", default=29674, help="The maximum position for site need to be masked.\nDefault: 29674")

    args = parser.parse_args()

    # Normolize the lineage features' mutaion
    dirpath = args.dirpath
    del_lin = args.delete_lineage
    len_UAB = int(args.threshold)
    geno_min = int(args.masked_min)
    geno_max = int(args.masked_max)

    fv_path = os.path.join(dirpath, '1_filtered_data', 'fv.txt')

    with open(dirpath + fv_path, "r") as f:
        for i in csv.reader(f):
            if i[0].startswith("X") or (i[2] == "") or (i[0] in del_lin):
                continue
            else:
                norm_mutate = normalize_position(i[2:], geno_min, geno_max)
                with open(dirpath + "1_filtered_data/fv_norm.txt", "a+") as fn:
                    fn.write(i[0] + "," + i[1] + "," + "".join(norm_mutate[1:-1].replace("'", "").split()) + "\n")

    # Normolize the samples' mutaion
    snps_path = "1_filtered_data/snp.txt"
    snp_num = 0
    snp_qualified = 0

    if (geno_min > 1) and (geno_max < 29903):
        with open(dirpath + snps_path, "r") as f:
            snp_num += 1
            for i in f.readlines():
                i_2 = i.strip().split(":")
                i_v = i_2[1].split(",")
                nor_snp = normalize_position(i_v, geno_min, geno_max)
                if nor_snp.count(",") + 1 < (2 * len_UAB):
                    continue
                else:
                    snp_qualified += 1
                    with open(dirpath + "1_filtered_data/snp_norm.txt", "a+") as fsn:
                        fsn.write(i_2[0] + ":" + "".join(nor_snp[1:-1].replace("'", "").split()) + "\n")

    elif (geno_min <= 1) and (geno_max >= 29903):
        with open(dirpath + snps_path, "r") as f:
            snp_num += 1
            for i in f.readlines():
                i_2 = i.strip().split(":")
                i_v = i_2[1].split(",")
                if i_v.count(",") + 1 < (2 * len_UAB):
                    continue
                else:
                    snp_qualified += 1
                    with open(dirpath + "1_filtered_data/snp_norm.txt", "a+") as fsn:
                        fsn.write(i_2[0] + ":" + "".join(i_v[1:-1].replace("'", "").split()) + "\n")

    print("Number of the qualified sequences: ", snp_num)
    print("------------------- Mutaions_normalized_completed ---------------------------")


if __name__ == '__main__':
    main()
