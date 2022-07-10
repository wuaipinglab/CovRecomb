'''
Function: Date acquisition and filteration, get mutations for each qualified sample and the lineage-defining feature mutations.

Input: Two files: the SARS-CoV-2 genomes and the corresponding epidemiology information from GISAID database, including the metadata.tsv and the mas.fasta.

Output: Two files: mutations for qualified samples (snp.txt) and the lineage-defining feature mutation (fv.txt).
'''

import re
import warnings
from multiprocessing import Pool
import pandas as pd
from Bio import SeqIO
import argparse
import os


def read_seq(file, meta, REF, times):
    print('Times: ', times)
    info_seqs = meta.index.tolist()

    msa_seqs = []
    seqs = {}
    with open(file) as f:
        for n, seq_record in enumerate(SeqIO.parse(f, 'fasta')):
            if seq_record.description != 'reference_sequence':
                seq_id = seq_record.description.split('|')[1]
                nstart = 500000 * (times - 1)
                nend = 500000 * times
                if nstart <= n < nend or seq_id == REF:
                    msa_seqs.append(seq_id)
                    seqs[seq_id] = str(seq_record.seq)
    print('Number of genomes: ', len(seqs))

    if nend > n:
        finish_all_seq = True
    else:
        finish_all_seq = False

    miss_info_seqs = list(set(msa_seqs) - set(info_seqs))
    for miss_info_seq in miss_info_seqs:
        seqs.pop(miss_info_seq)
    print('Number of genomes with info: ', len(seqs))

    ambiguous_time_seqs = []
    for seq_id in seqs:
        if not re.search('20\d\d-\d\d-\d\d', meta.loc[seq_id, 'date']):
            ambiguous_time_seqs.append(seq_id)
    for ambiguous_time_seq in ambiguous_time_seqs:
        seqs.pop(ambiguous_time_seq)
    print('Number of genomes with complete info: ', len(seqs))

    return finish_all_seq, seqs


def get_genome_sites(sequence):
    sites = []
    n = 0
    for i in sequence:
        if i != '-':
            sites.append(n)
        n += 1
    return sites


def trim_seq(sequence, sites):
    seq_trimed_list = []
    for i in sites:
        seq_trimed_list.append(sequence[i])
    seq_trimed = ''.join(seq_trimed_list)
    return seq_trimed


def is_qualified(sequence, length):
    if sequence.count('A') + sequence.count('T') + sequence.count('G') + sequence.count('C') >= length:
        return True
    else:
        return False


def clean_seq(seq_id, sequence, sites, ref_genome, length):
    seq_trimed = trim_seq(sequence, sites)
    if is_qualified(seq_trimed, length):
        snps = get_snps(seq_trimed, ref_genome)
        return {seq_id: snps}
    return False


def get_snps(sequence, ref_genome):
    snps = []
    gaps = []
    for gsite in range(0, len(ref_genome)):
        if sequence[gsite] != ref_genome[gsite] and sequence[gsite] in ['A', 'T', 'G', 'C']:
            snps.append(str(gsite + 1) + '_' + sequence[gsite])
        elif sequence[gsite] == '-':
            gaps.append(gsite + 1)

    continuous_gaps = []
    for gap in gaps:
        if len(continuous_gaps) == 0 or gap == continuous_gaps[-1] + 1:
            continuous_gaps.append(gap)
        else:
            if len(continuous_gaps) % 3 == 0 and continuous_gaps[0] != 1 and continuous_gaps[-1] != 29891:
                snps.append(str(continuous_gaps[0]) + '_' + '-' * len(continuous_gaps))
            continuous_gaps = [gap]

    return snps


def write_new_file(file, seq_snps):
    with open(file, 'w') as f:
        for i in seq_snps:
            f.write(i + ':' + ','.join(seq_snps[i]) + '\n')


def creat_dir(new_dir):
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)


def main():
    # command line interface
    parser = argparse.ArgumentParser(description='Read the qulified sequences.', usage='''python3 fasta_clean_get_FM.py''')

    parser.add_argument("-dir", "--dirpath", default="/home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12/", help="The totoal folder address of this work.\nDefault: /home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12")
    parser.add_argument("-m", "--meta_path", default='0_raw_data/metadata-2.tsv', help="The subfolder address of metadata.tsv.\nDefault: 0_raw_data/metadata-2.tsv")
    parser.add_argument("-s", "--sequence_path", default='0_raw_data/msa_0220/msa_0220.fasta', help="The subfolder address of the sequences.fasta.\nDefault: 0_raw_data/msa_0220/msa_0220.fasta")
    parser.add_argument("-r", "--reference", default="EPI_ISL_402125", help="The reference sequence ID, indicating which sequence\nin the alignment is the reference (by sequence ID).\nDefault: EPI_ISL_402125")
    parser.add_argument("-c", "--core_number", default=8, help="The number of cores used while computation running.\nDefault: 8")
    parser.add_argument("-l", "--sequence_length", default=27000, help="The least number of the 'ATCG' for a qualified sequence.\nDefault: 27000")
    parser.add_argument("-fm", "--feature_mutation", default=0.75, help="The percentage of shared mutations between a lineage could be considered as the feature mutations.\nDefault: 0.75")
    args = parser.parse_args()

    dirpath = args.dirpath
    REF = args.reference
    meta_path = dirpath + args.meta_path
    sequence_path = dirpath + args.sequence_path
    cor_num = int(args.core_number)
    length = int(args.sequence_length)
    fm_threshold = float(args.feature_mutation)

    creat_dir(dirpath + '1_filtered_data')
    snps_path = dirpath + '1_filtered_data/snp.txt'

    warnings.filterwarnings('ignore')

    meta = pd.read_csv(meta_path, delimiter='\t', index_col=2).dropna(axis=0, subset=['date', 'region_exposure', 'country_exposure', 'division_exposure'])
    meta = meta[(meta['host'] == 'Human')]

    times = 0
    finish_all_seq = False
    seq_snps = {}
    while not finish_all_seq:
        times += 1
        finish_all_seq, seqs = read_seq(sequence_path, meta, REF, times)
        sites = get_genome_sites(seqs[REF])
        ref_genome = trim_seq(seqs[REF], sites)
        result_list = []
        p = Pool(cor_num)
        for i in seqs:
            result = p.apply_async(clean_seq, args=(i, seqs[i], sites, ref_genome, length))
            result_list.append(result)
        p.close()
        p.join()
        for result in result_list:
            if result.get():
                seq_snps.update(result.get())

    seq_snps = dict(sorted(seq_snps.items(), key=lambda x: int(x[0].split('_')[2])))
    print('Number of genomes with high quality: ', len(seq_snps))

    write_new_file(snps_path, seq_snps)

    fv_path = os.path.join(dirpath, '1_filtered_data', 'fv.txt')

    # PANGO_WHO = {
    #     'B.1.1.7': 'Alpha',
    #     'Q': 'Alpha',
    #     'B.1.351': 'Beta',
    #     'P.1': 'Gamma',
    #     'B.1.617.2': 'Delta',
    #     'AY': 'Delta',
    #     'C.37': 'Lambda',
    #     'B.1.621': 'Mu'
    # }

    meta = pd.read_csv(meta_path, delimiter='\t', index_col=2).dropna(axis=0, subset=['pango_lineage'])

    lineages = {}
    with open(snps_path) as f:
        lines = f.readlines()
        for line in lines:
            seq_id = line.split(':')[0]
            if seq_id in meta.index:
                lineage = meta.loc[seq_id, 'pango_lineage']

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
            # if lineages[lineage]['count'] / len(meta) > 0.0001:
            fv[lineage] = []
            for snp in lineages[lineage]['snps']:
                if lineages[lineage]['snps'][snp] / lineages[lineage]['count'] >= fm_threshold:
                    fv[lineage].append(snp)

    with open(fv_path, 'w') as f:
        for lineage in fv:
            if '' in fv[lineage]:
                fv[lineage].remove('')
            fv[lineage] = sorted(fv[lineage], key=lambda x: int(x.split('_')[0]))
            f.write(lineage + ',' + str(lineages[lineage]['count']) + ',' + ','.join(fv[lineage]) + '\n')

    print("------------------- Sequences_filteration_and_mutations_calling_completed ---------------------------")


if __name__ == '__main__':
    main()
