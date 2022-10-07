import csv
from collections import Counter

def get_keys(d, value):
    return [k for k,v in d.items() if v == value]


def merge_dict(x,y):
    for k,v in x.items():
        if k in y.keys():
            y[k] += v
        else:
            y[k] = v


DIRPATH = infile = "/home/soniali/Desktop/03_CovRecomb/CovRecomb/MERS/"
lineage_file = 'fv_75.txt'
Lineage_v = {}
linA_list = []
feature_mutations = []
with open(DIRPATH+lineage_file,'r') as f:
    count = 0
    for i in csv.reader(f):
        if i[2] != "":
            if len(i[2:])>= 1:
                linA_list.append(i[0])
                Lineage_v[i[0]] = i[2:]

Strain_list_snp = []
variants_all = {}
with open(DIRPATH+"snp.txt","r") as f:
    for i in f.readlines(): 
        i_2 = i.split(':')  
        Strain_list_snp.append(i_2[0])
        variants_all[i_2[0]] = i_2[1].strip().split(',')

Strain_list_meta = []
pango_lineage = {}
with open(infile+"group.csv","r") as f:
    next(f)
    for i in csv.reader(f,delimiter=','):
        Strain_list_meta.append(i[2])
        pango_lineage[i[2]] = i[1]

raw_count={'0':0,'1':0}
for lin in linA_list:
    lin_denovo_num = []

    lin_epi = [k for k,v in pango_lineage.items() if v == lin]
    for i in lin_epi:
        epi_record = []
        num_denovo = 0
        try:
            for v in variants_all[i]:
                if v in set(Lineage_v[lin]):
                    epi_record.append(num_denovo)
                    num_denovo = 0
                else:
                    num_denovo += 1
            epi_record.append(num_denovo)
            denovo_count = dict(Counter(epi_record))
            raw_count = dict(Counter(raw_count)+Counter(denovo_count))
        except KeyError:
            continue

lin_denovo_count = dict([(k,raw_count[k]) for k in sorted(raw_count.keys())])
num = 0
for d in lin_denovo_count:
    if d != 0:
        num = num + lin_denovo_count[d]
    else:
        continue

lin_denovo_prop = {}
for d in lin_denovo_count:
    if d != 0 :
        lin_denovo_prop[d] = lin_denovo_count[d]/num
    else:
        continue

Proportion = list(lin_denovo_prop.keys())[0:13]
CF = list(lin_denovo_prop.values())[0:13]
# [1, 2, 3, 4]
# [0.7108108108108108, 0.21621621621621623, 0.05675675675675676, 0.016216216216216217]
