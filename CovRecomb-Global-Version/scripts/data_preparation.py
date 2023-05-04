import re
import warnings
import os

print("\n", "---------- Step1: Get the mutations and deletions for each qualified sample  ----------", "\n")
###  import the essential variables and parameters
dirpath = "/home/soniali/Desktop/02_recom_230203/data/2023_02_08/"
os.chdir(dirpath+"scripts/")
from function_set import snps_path, fm_threshold, length, geno_min, geno_max, dirpath, meta_path, sequence_path, parameter_souce, fv_feature_filename, \
    sort_humanly, creat_dir, write_new_file, cluster_similar_lineages, min_continous_mutations

creat_dir(dirpath+'/1_filtered_data')
warnings.filterwarnings('ignore')

###  read the meta information for each genome
with open(meta_path, "r") as f:
    next(f)
    rows = f.readlines()

pango_lineage = {}
collect_date = {}
for i in rows:
    info = i.split("\t")
    seq_length = int(info[8])
    host = info[9]
    if host == "Human" and seq_length >= length and re.search('20\d\d-\d\d-\d\d', info[5]):
        epi = info[0]
        pango_lineage[epi] = info[13]
        collect_date[epi] = info[5]
        
del rows
print("Number of qualified genomes in meta file: {}".format(len(pango_lineage))) # 14206354 ////14810273 / 7695697

###  read the mutations information for each genome
with open(sequence_path, "r") as h:
    next(h)
    hrows = h.readlines()

variants_all = {}
for i in hrows:
    info = i.split("\t")
    if len(info[19]) > 0:
        epi = info[0].split("|")[0]
        mut = info[19].split(",")
        dele = info[20].split(",")
        mut_list = []
        for v in mut:
            if geno_min <= int(v[1:-1]) <= geno_max:
                mut_list.append(v[1:-1]+"_"+v[-1])
        del_list = []
        for v in dele:
            try:
                startpos = v.split("-")[0]
                if geno_min <= int(startpos) <= geno_max:
                    endpos = v.split("-")[1]
                    del_len = startpos+"_"+"-"*(int(endpos) - int(startpos) +1)
                    del_list.append(del_len)
            except:
                continue
        all_mut = sort_humanly(del_list + mut_list)
        variants_all[epi] = all_mut

del hrows
write_new_file(snps_path, variants_all)
print('Number of qualified genomes in nextclade file: ', len(variants_all))

### ---------------------------Construction-------------------------------
print("\n", "---------- Step2: Get the feature mutations for each lineage ----------", "\n")
###  calculate the feature mutations for each lineage
fv_path = os.path.join(dirpath, '1_filtered_data','fv.txt')
lineages = {}
with open(snps_path) as f:
    lines = f.readlines()
    for line in lines:
        seq_id = line.split(':')[0]
        if seq_id in pango_lineage:
            lineage = pango_lineage[seq_id]
            
            if lineage not in lineages:
                lineages[lineage] = {'count': 1, 'snps': {}}
            else:
                lineages[lineage]['count'] += 1
            snps = line.strip().split(':')[1].split(',')
            for snp in snps:
                lineages[lineage]['snps'][snp] = lineages[lineage]['snps'].get(snp, 0) + 1

lineages = dict(sorted(lineages.items(), key=lambda x: x[0]))

fv = {}
for lineage in lineages:
    if lineage != 'None':
        fv[lineage] = []
        for snp in lineages[lineage]['snps']:
            if lineages[lineage]['snps'][snp] / lineages[lineage]['count'] >= fm_threshold:
                fv[lineage].append(snp)

with open(fv_path, 'w') as f:
    for lineage in fv:
        if '' in fv[lineage]:
            fv[lineage].remove('')
        fv[lineage] = sorted(fv[lineage], key=lambda x: int(x.split('_')[0]))
        f.write(lineage+','+str(lineages[lineage]['count'])+','+','.join(fv[lineage])+'\n')

print("\n", "---------- Step2: Calculate the parameter for continuous feature mutaiton and cluster similar lineages --------------", "\n")

if parameter_souce == "Self-define":
    import csv
    len_UXY = min_continous_mutations(fv_path, pango_lineage, variants_all)
    Lineage_v = {}
    lineage_n = {}
    with open(fv_path, 'r') as f:
        for i in csv.reader(f):
            if (i[0].startswith("X")) or (i[0] == "Unassigned") or (i[0] == "None") or (len(i[2:]) < len_UXY) or (int(i[1]) <= 100):
                continue
            else:
                lineage_n[i[0]] = i[1]
                Lineage_v[i[0]] = i[2:]
    with open(dirpath + "1_filtered_data/" + fv_feature_filename, "a+") as h:
        for lin in lineage_n:
            h.write(lin + "," + str(lineage_n[lin]) + "," + ",".join(Lineage_v[lin]) + "\n")

elif parameter_souce == "Default":
    from function_set import len_UXY
    mini_simi = cluster_similar_lineages(fv_path, len_UXY, fv_feature_filename)
    with open(dirpath+"scripts/function_set.py","a+") as f:
        f.write("\n"+"mini_simi_default = "+ str(mini_simi)+"\n")

print("------------------- Preparation for qualified genomes and lineages are completed ---------------------------")