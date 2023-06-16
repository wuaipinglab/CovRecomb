import os
import subprocess
import datetime
import pandas as pd

def TimeConverter(ms):
    s = round(ms / 1000, 2)
    m = round(s / 60, 2)
    h = round(m / 60, 2)
    return s, m, h


cor_num = 2

dirpath = "/HOME/usher/data/jiayingli/"
for n in [0,1,2,3,4]:
    turns_path = dirpath+"2_first_test/turns"+str(n+1)+"/"
    os.chdir(turns_path+"rivet/")
    for file in os.listdir(turns_path):
        if file.startswith("real_recom_lineage_samples_"):
            child_fasta_file = turns_path+file
            
    # read sequences
    seq_samples = {}
    with open(child_fasta_file) as h:
        for s in h.readlines():
            if s[0] == ">":
                seq_id = s.strip().split(">")[1]
            else:
                seq_samples[seq_id] = str(s.strip())
                
    recom_total_num = sample_size = len(seq_samples)
    
    recom_path = turns_path+"rivet/real_recom_lineage_samples_rivet.fasta"
    with open(recom_path,"w+") as f:
        with open(child_fasta_file) as h:
            for s in h.readlines():
                if s[0] == ">":
                    seq_id = s.strip().replace(" ","_")
                    f.write(seq_id+"\n")
                else:
                    f.write(str(s.strip())+"\n")

    ref_path = dirpath+"0_raw_data/wuhCor1.fa"
    sars_path = dirpath+"0_raw_data/problematic_sites_sarsCov2.vcf"
    public_path = dirpath+"0_raw_data/public-latest.all.masked.pb.gz"

    print("rivet Start")
    starttime = datetime.datetime.now()
    subprocess.call (["cd %s && \
        echo y | mafft --thread 2 --auto --keeplength --addfragments %s %s > recom_aligned_seqs.fa" % (turns_path+"rivet/",recom_path,ref_path)],shell=True)

    subprocess.call (["cd %s && \
        echo y | faToVcf -maskSites=%s recom_aligned_seqs.fa recom_aligned_seqs.vcf" % (turns_path+"rivet/",sars_path)],shell=True)

    subprocess.call (["cd %s && \
        echo y | usher -T 2 -i %s -v recom_aligned_seqs.vcf -o recom_user_seqs.pb" % (turns_path+"rivet/",public_path)],shell=True)

    subprocess.call (["cd %s && \
        echo y | mkdir USER_SAMPLES_rivet" % (turns_path+"rivet/")],shell=True)

    subprocess.call (["cd %s && \
        echo y | grep -e '>' %s | perl -pi -e 's/>//' > recom_USER_SAMPLES_rivet.txt" % (turns_path+"rivet/",recom_path)],shell=True)
    
    subprocess.call (["cd %s && \
        echo y | ripples-fast -i recom_user_seqs.pb -s recom_USER_SAMPLES_rivet.txt -d USER_SAMPLES_rivet/ -T 2" % (turns_path+"rivet/")],shell=True)
    
    endtime = datetime.datetime.now()
    time_microsecond2 = (endtime -starttime).seconds * 1000 + (endtime -starttime).microseconds / 1000
    time_second2, time_minute2, time_hour2 = TimeConverter(time_microsecond2)
    
    print("rivet Results")
    import os
    ## record the descendants for each node
    with open(turns_path+"rivet/USER_SAMPLES_rivet/"+"descendants.tsv") as f:
        next(f)
        rows = f.readlines()
        
    node_des = {}
    for row in rows:
        node_des[row.strip().split("\t")[0]] = row.strip().split("\t")[1].split(",")

    ## find the recombinant node
    with open(turns_path+"rivet/USER_SAMPLES_rivet/"+"recombination.tsv") as f:
        next(f)
        rows = f.readlines()

    recom_node = []
    for row in rows:
        score1,score2 = row.strip().split("\t")[10], row.strip().split("\t")[11]
        if int(score1)>=3 and int(score2) >=3:
            recom_node.append(row.strip().split("\t")[0])

    ## record the analyzing samples
    seq_id = []
    with open(recom_path) as f:
        for i in f.readlines():
            if i.startswith(">"):
                seq_id.append(i.strip().split(">")[1])

    ## detect the recombinant samples    
    recom_id = []
    for node in recom_node:
        recom_id.extend(list(set(node_des[node]) & set(seq_id)))

    ## calculat the tpr for rivet
    final_recom_id = set(recom_id)
    TPR = round(len(final_recom_id)/len(seq_id),4)
    
    rivet_results = [str(recom_total_num),str(TPR), str(time_microsecond2), str(time_second2), str(time_minute2), str(time_hour2),str(cor_num)]
    print(rivet_results)

    # Write results to output file
    with open(turns_path+"methods_results_0615.csv","a+")as h:
        h.write("turns"+str(n+1)+","+"rivet"+","+",".join(rivet_results)+"\n")

    print("rivet Finished")

    
# nohup python3 ./RIPPLES_run.py > RIPPLES_run.log 2>&1 &

## Examples:
# mafft --thread 10 --auto --keeplength --addfragments test_samples.fa wuhCor1.fa > aligned_seqs.fa
# faToVcf -maskSites=problematic_sites_sarsCov2.vcf aligned_seqs.fa aligned_seqs.vcf
# usher -T 10 -i public-latest.all.masked.pb.gz -v aligned_seqs.vcf -o user_seqs.pb
# mkdir USER_SAMPLES_rivet/
# grep -e '>' test_samples.fa | perl -pi -e 's/>//' > USER_SAMPLES_rivet.txt
# rivet -i user_seqs.pb -s USER_SAMPLES_rivet.txt -d USER_SAMPLES_rivet/ -T 60
