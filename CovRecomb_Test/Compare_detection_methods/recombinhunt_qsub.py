import csv
import datetime

import json
import pandas as pd

from recombinhunt_exp import cal, compare_lists, recombinhunt_experiment, time_converter

# 运行recombinhunt
performance_turns_results = []
performance_turns_times = []
speed_sample_size = [72, 144, 216, 360]
for i in range(1, 7):
    start_time = datetime.datetime.now()
    results_dict = recombinhunt_experiment(f'performance/turns{i}/real_recom_lineage_samples_72.fasta')
    # results_dict = recombinhunt_experiment(f'speed/real_recom_lineage_samples_{speed_sample_size[i-1]}.fasta')
    end_time = datetime.datetime.now()
    delta = end_time - start_time
    performance_turns_results.append(results_dict)
    performance_turns_times.append(delta)

# 计算覆盖率和准确率所需其他数据
with open('../data/raw/alias_key.json', 'r', encoding='utf-8') as file:
    recombinant_lineages = json.load(file)

# 计算性能指标并写入csv文件
with open('../results/test.csv','w') as csvfile:
    writer = csv.writer(csvfile)

    writer.writerow(['turns','method','sample_size','tpr','time_microsecond',
                     'time_second','time_minute','time_hour','accuracy_of_parental_lineages', 'core'])
    for i in range(1, 7):
    # for i in range(1, 5):
        meta = pd.read_csv(f'../data/raw/performance/turns{i}/real_recom_lineage_samples_72.fasta',
        # meta = pd.read_csv(f'../data/intermediate/speed/real_recom_lineage_samples_{speed_sample_size[i-1]}.fasta.tsv',
                           sep='\t', usecols=['seqName', 'Nextclade_pango'])
        coverage, accuracy = cal(performance_turns_results[i-1], recombinant_lineages, meta)
        time_microsecond = performance_turns_times[i-1].seconds * 1000 + performance_turns_times[i-1].microseconds / 1000
        time_second, time_minute, time_hour = time_converter(time_microsecond)
        results= [f'turns{i}', 'recombinhunt', str(len(meta)), str(coverage), str(time_microsecond),
                  str(time_second), str(time_minute), str(time_hour), str(accuracy), '1']
        writer.writerow(results)
        