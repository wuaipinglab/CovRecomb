from tqdm import tqdm

import pandas as pd
from recombinhunt.core.environment import Environment
from recombinhunt.core.method import Experiment
from recombinhunt.validation.utils import LineageHierarchy
import subprocess


def nextclade_2_recombinhunt(substitutions, deletions):
    """
    Converts nextclade's output to recombinhunt's input format
    """
    substitutions = [] if pd.isna(substitutions) else substitutions.split(',')
    deletions = [] if pd.isna(deletions) else deletions.split(',')
    mutation_list = []
    deletion_list = []
    for mutation in substitutions:

        if mutation:
            pos = ''.join(filter(str.isdigit, mutation))
            ref = mutation[0]
            alt = mutation[-1]
            mutation_list.append(f'{pos}_{ref}|{alt}')

    for deletion in deletions:
        if deletion:
            deletion = deletion.replace('-', '_')
            deletion_list.append(deletion)

    combined_list = mutation_list + deletion_list
    combined_list.sort(key=lambda x: int(x.split('_')[0]))
    return combined_list


def recombinhunt_experiment(fasta_file):
    """
    Runs recombinhunt on a given fasta file
    """
    env = Environment('../models/env_nextstrain_2023_03_30')
    lh = LineageHierarchy("../models/alias_key.json")
    custom_environments = {
        "XBB": Environment('../models/env_nextstrain_2023_03_30', ignore_lineage=[l for l in env.included_lineages() if
                                                                                  l.startswith("XBB") or any(
                                                                                      [sl in l for sl in (
                                                                                          "EG", "EK", "EL", "EM", "EU",
                                                                                          "FD", "FE", "FG", "FH",
                                                                                          "FL")])]),
        "XBF": Environment('../models/env_nextstrain_2023_03_30',
                           ignore_lineage=[l for l in env.included_lineages() if l.startswith("XBF")]),
        "XAY": Environment('../models/env_nextstrain_2023_03_30',
                           ignore_lineage=[l for l in env.included_lineages() if l.startswith("XAY")]),
        "XP": Environment('../models/env_nextstrain_2023_03_30',
                          ignore_lineage=[l for l in env.included_lineages() if l.startswith("XBB")])
    }

    _ = subprocess.run(
        f'nextclade run --input-dataset ../models/nextclade --output-tsv=../data/intermediate/{fasta_file}.tsv ../data/raw/{fasta_file}',
        shell=True)

    nextclade_df = pd.read_csv(f'../data/intermediate/{fasta_file}.tsv', sep='\t',
                               usecols=['seqName', 'substitutions', 'deletions', 'Nextclade_pango'])

    results_dict = dict()
    for _, row in tqdm(nextclade_df.iterrows(), total=len(nextclade_df)):
        recombinhunt_list = nextclade_2_recombinhunt(row['substitutions'], row['deletions'])

        exp = Experiment(custom_environments.get(row['Nextclade_pango'], env), lh)
        exp.set_target(recombinhunt_list)
        try:
            result = exp.run()
        except Exception as e:
            result = None
            print(f'Error in {row["seqName"]}: {e}')
        results_dict[row['seqName']] = result
    return results_dict


def compare_lists(pattern_list, obj_list):
    if len(pattern_list) != len(obj_list):
        return False

    def match_with_wildcard(pattern, obj):
        if pattern[-1] == '*':
            prefix = pattern[:-1]
            if not obj.startswith(prefix):
                return False
            # 检查剩余部分是否只包含版本号
            remaining = obj[len(prefix):]
            return remaining == '' or (remaining.startswith('.') and all(part.isdigit() for part in remaining.split('.')[1:]))
        return pattern == obj

    def backtrack(pattern_index, used_objects):
        if pattern_index == len(pattern_list):
            return len(used_objects) == len(obj_list)

        pattern = pattern_list[pattern_index]
        for i, obj in enumerate(obj_list):
            if i not in used_objects and match_with_wildcard(pattern, obj):
                used_objects.add(i)
                if backtrack(pattern_index + 1, used_objects):
                    return True
                used_objects.remove(i)

        return False

    return backtrack(0, set())


def cal(result, recombinant_lineages, meta):
    tp = 0
    right_in_dataset = 0
    for key, val in result.items():
        if key is None:
            continue
        if len(val.p_values) > 1:
            tp += 1
            nextclade_pango = meta.loc[meta['seqName'] == key, 'Nextclade_pango'].values[0].split('.')[0]

            ture_parents = recombinant_lineages[nextclade_pango]
            parents = list(set([_[1] for _ in val.p_values.keys()]))
            if compare_lists(list(set(ture_parents)), parents):
                right_in_dataset += 1
    coverage = tp / len(result)

    if coverage == 0:  # 如果coverage是0则无需计算accuracy
        return 0, None

    accuracy = right_in_dataset / tp

    return coverage, accuracy


def time_converter(ms):
    # 保留两位小数，但若ms太小，h就会显示为0。
    s = round(ms / 1000, 2)
    m = round(s / 60, 2)
    h = round(m / 60, 2)
    return s, m, h
