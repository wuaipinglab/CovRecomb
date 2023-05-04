import pandas as pd
import copy
import os
import csv
import re
import datetime
import numpy as np

def average(arg):
    return round(sum(arg)/len(arg), 4)


def get_linXY_mean_std(linX_name,lineage,Y_muts,linY_name,X_muts):
    
    mutY_diff_prop = verify_linXY_prop_mean_std(linX_name,lineage,Y_muts,linY_name)
    mutX_diff_prop = verify_linXY_prop_mean_std(linY_name,lineage,X_muts,linX_name)
    Confidence_mean = round((average(mutY_diff_prop)+average(mutX_diff_prop))/2,4)
    
    return Confidence_mean
    
def verify_linXY_prop_mean_std(linX_name,lineage,Y_muts,linY_name):

    mutY_linX_prop = []
    for snp in Y_muts:
        try:
            mutY_linX_prop.append(round(lineage[linY_name]['snps'][snp] / lineage[linY_name]['count'],4))
        except:
            mutY_linX_prop.append(0)

    mutX_linX_prop = []
    for snp in Y_muts:
        try:
            mutX_linX_prop.append(round(lineage[linX_name]['snps'][snp] / lineage[linX_name]['count'],4))
        except:
            mutX_linX_prop.append(0)
    
    mutY_diff_prop = [round(mutY_linX_prop[i] - mutX_linX_prop[i],4) for i in range(len(Y_muts))]
    
    return mutY_diff_prop


def main():
    ###  import the essential parameters
    dirpath = "/home/soniali/Desktop/02_recom_230203/data/2023_02_08/"
    os.chdir(dirpath+"scripts/")
    from function_set import fv_feature_filename, snps_path, meta_path, length, mini_simi_default, parameter_souce, \
        creat_dir,del_star, input_meta_info, detect_mode,diff_days_fm,get_keys,extract_AB_bk,find_AB_ances, \
            get_time_epi,left_right_lin,find_bk, summary_epi_dist, sort_humanly, calcul_bk,calcul_bk_2

    ############## part1  read snp for each sample
    print("\n","----------------------- Load datasets -----------------------","\n")
    lineage_file = dirpath + "1_filtered_data/" + fv_feature_filename
    dirpath_pattern = dirpath + "3_recom_pattern/"
    creat_dir(dirpath_pattern)
    dirpath_recomb = dirpath + "2_recomb_identified/"
    verified_output = dirpath_recomb+"0_putative_recombinants_step2.csv"
    path_breakloc = dirpath_recomb + "bk_info/"
    
    ## read each sample's mutations
    variants_all = {}
    with open(snps_path, "r") as f:
        for i in f.readlines():
            i_2 = i.split(':')
            variants_all[i_2[0]] = i_2[1].strip().split(',')
            
    ## read the meta file
    with open(meta_path, "r") as f:
        next(f)
        rows = f.readlines()

    num_meta = 0
    pango_lineage = {}
    collect_date = {}
    country = {}
    region = {}
    for i in rows:
        num_meta += 1
        info = i.split("\t")
        seq_length = int(info[8])
        host = info[9]
        if host == "Human" and seq_length >= length and re.search('20\d\d-\d\d-\d\d', info[5]):
            epi = info[0]
            region[epi] = info[6].split("/")[0].strip()
            country[epi] = info[6].split("/")[1].strip()
            pango_lineage[epi] = info[13]
            collect_date[epi] = info[5]

    del rows
    if parameter_souce == "Default":
        mini_simi = mini_simi_default
        mini_simi = mini_simi_default
        cluster_dict = {}
        with open(dirpath+"/1_filtered_data/cluster"+str(mini_simi)+".csv","r") as f:
            cluster = csv.reader(f)
            for i in cluster:
                cluster_dict[i[0]] = i[1:-1]

        Lineage_v = {}
        with open(lineage_file,'r') as f:
            for i in csv.reader(f):
                if "*" not in i[0]:
                    Lineage_v[i[0]] = i[2:]
                else:
                    lin_name = i[0].split("*")[0]
                    for l in cluster_dict:
                        if lin_name in cluster_dict[l]:
                            cluster_n = l
                    for c in cluster_dict[cluster_n]:
                        Lineage_v[c] = i[2:]

    elif parameter_souce == "Self-define":
        Lineage_v = {}
        with open(lineage_file, 'r') as f:
            for i in csv.reader(f):
                Lineage_v[i[0]] = i[2:]

    print("\n","----------------------- Find recombination events -----------------------","\n")
    df = pd.read_csv(verified_output)
    for epi in list(df["sample_id"]):
        find_bk(df, epi, variants_all, path_breakloc)

    ## Distinguish ancient-like recombination events to the transmitted genomes
    df= df.sort_values(by = ["sample_id","mutation_pattern"],ascending = ("False","False"))
    df["situation"] = df['lineage_X']+"_"+ df['lineage_Y']
    df_final = df.sort_values(by = ["collect_date"],ascending = ("True"))
    df_grouped = pd.DataFrame(columns = list(df_final))
    each_group = {}
    for i in set(df_final['situation']): 
        temp_ABC = (df_final.loc[df_final["situation"] == i]).sort_values(by = ["collect_date"],ascending = ("True"))
        df_first = pd.DataFrame(columns = list(df_final))
        df_candi = pd.DataFrame(columns = list(df_final))
        
        if temp_ABC.shape[0] == 1:
            df_first = pd.concat([df_first,temp_ABC])
            df_grouped = pd.concat([df_grouped,df_first])
            group_epi = list(temp_ABC["sample_id"])[0]
            each_group[group_epi] = [group_epi]
            
        else:
            df_first = pd.concat([df_first,temp_ABC.head(1)])
            ABC_index = (temp_ABC.index).tolist()
            left_A = del_star(temp_ABC.iloc[0]["lineage_X"])
            right_B = del_star(temp_ABC.iloc[0]["lineage_Y"])
            epi = temp_ABC.iloc[0]["sample_id"]

            first_mode = detect_mode(left_A, right_B, variants_all, epi, Lineage_v)
            turn_mode = [first_mode]
            for n in ABC_index[1:]:
                temp_mode = detect_mode(left_A,right_B,variants_all,temp_ABC.loc[n]["sample_id"], Lineage_v)
                turn_mode2 = copy.deepcopy(turn_mode)
                former_time = temp_ABC.loc[ABC_index[ABC_index.index(n)-1]]["collect_date"]
                temp_time = temp_ABC.loc[n]["collect_date"]
                if (temp_time <= (datetime.datetime.strptime(former_time, "%Y-%m-%d") + datetime.timedelta(days=+30)).strftime("%Y-%m-%d")):
                    turn_mode,df_candi = diff_days_fm(turn_mode2,df_candi,temp_mode,temp_ABC,n,turn_mode,4)
                else:
                    df_candi = df_candi.append(temp_ABC.loc[n])
                    turn_mode.append(temp_mode)

            ## group offspring to the identified events
            df_situation = pd.concat([df_first,df_candi])
            df_grouped = pd.concat([df_grouped,df_situation])
            
            all_epi = list(temp_ABC.sample_id)
            if df_situation.shape[0] == 1:
                group_epi = list(df_first["sample_id"])[0]
                each_group[group_epi] = all_epi
            else:
                recom_event = list(df_situation["sample_id"])
                recom_event_mut = {}
                for ev in recom_event:
                    recom_event_mut[ev] = variants_all[ev]
                
                already_ev = []
                for temp_id in all_epi:
                    if temp_id in recom_event:
                        already_ev.append(temp_id)
                        each_group[temp_id] = [temp_id]
                    else:
                        epi_mut = variants_all[temp_id]
                        inher_numt_num = {}
                        for cand in already_ev:
                            inher_numt_num[cand] = len(set(recom_event_mut[cand]) & set(epi_mut))
                        
                        res = []
                        for key, value in inher_numt_num.items():
                            if value == max(inher_numt_num.values()):
                                res.append(key)
                        
                        if len(res) == 1:
                            each_group[res[0]].extend([temp_id])
                        else:
                            ca_num = {}
                            for ca in res:
                                ca_num[ca] = len(set([region[ca],country[ca]]) & set([region[temp_id],country[temp_id]]))
                        
                            cas = []
                            for key, value in ca_num.items():
                                if value == max(ca_num.values()):
                                    cas.append(key)
                            
                            if len(cas) == 1:
                                each_group[cas[0]].extend([temp_id])
                            else:
                                each_group[cas[-1]].extend([temp_id])

    each_group_sort = dict(sorted(each_group.items(), key=lambda item:len(item[1]),reverse=True))

    df_grouped_sort = df_grouped.sort_values(by =["collect_date"] ,ascending = ("True"))
    df_grouped_sort = input_meta_info(df_grouped_sort, region, country)
    df_grouped_sort.to_csv(dirpath_pattern+"Independent_recombination_event_"+str(df_grouped_sort.shape[0])+".csv",index=None)

    with open(dirpath_pattern+"each_group_sort_"+str(df_grouped_sort.shape[0])+".txt","a+") as f:
        events_num = 0
        for n in each_group_sort:
            if len(each_group_sort[n]) > 1:
                events_num += 1
            f.write(n+":"+str(len(each_group_sort[n]))+","+",".join(each_group_sort[n])+"\n")

    with open(dirpath_pattern+"each_group_sort_"+str(events_num)+".txt","a+") as f:
        for n in each_group_sort:
            if len(each_group_sort[n]) > 1:
                f.write(n+":"+str(len(each_group_sort[n]))+","+",".join(each_group_sort[n])+"\n")

    ## Delete those events with single genome
    df_final = pd.DataFrame(columns=df_grouped_sort.columns)
    delete_id = []
    for i in df_grouped_sort.index:
        n = df_grouped_sort.loc[i, "sample_id"]
        if len(each_group_sort[n]) == 1:
            delete_id.append(n) 
        else:
            df_final = df_final.append(df_grouped_sort.loc[i,])

    for i in df_final.index:
        n = df_final.loc[i, "sample_id"]
        num = int(len(each_group_sort[n]))
        
        epi_country = []
        for id in each_group_sort[n]:
            epi_country.append(country[id])
        stast_country = {}
        for c in epi_country:
            stast_country[c] = int(epi_country.count(c))

        geo_dis = ""
        for a in stast_country:
            geo_dis += a+"_"+str(stast_country[a])+"/"
            
        df_final.loc[i,"Number_of_epidemic_recombinant"] = round(num)
        df_final.loc[i,"Geographical_distribution"] = geo_dis
    df_final = df_final.sort_values(by = ["Number_of_epidemic_recombinant"],ascending = False)
    df_final = input_meta_info(df_final, region, country)
    df_final.to_csv(dirpath_pattern+"Independent_recombination_event_"+str(df_final.shape[0])+".csv",index=None)

    # df_final = pd.read_csv(dirpath_pattern+"Independent_recombination_event_1451.csv")
    within_events = df_final["sample_id"].tolist()
    with open(dirpath_pattern+"each_group_sort_4888.txt","r") as f:
        rows = f.readlines()
        
    time_delete_epi = []
    for i in rows:
        epi = i.split(":")[0]
        if set([epi]) - set(within_events) != set():
            time_delete_epi.extend(i.strip().split(":")[1].split(",")[1:])
            
    len(time_delete_epi)
    df_putative_final = pd.read_csv(dirpath + "2_recomb_identified/0_putative_recombinants_104004.csv")
    df_putative_final2 = df_putative_final.loc[~df_putative_final["sample_id"].isin(time_delete_epi)]
    df_putative_final2 = input_meta_info(df_putative_final2, region, country)
    df_putative_final2.to_csv(dirpath + "2_recomb_identified/1_putative_recombinants_final_"+str(df_putative_final2.shape[0])+".csv",index = None)
    
    ## Divide the identified recombination events into "X" series or not.
    # df = pd.read_csv(dirpath + "2_recomb_identified/1_putative_recombinants_final_135567.csv")
    df = df_putative_final2
    X_series = {}
    x_num = 0
    for i in df.index:
        pan = df.loc[i,"pango_lineage"]
        if pan[0] == "X" and pan not in X_series:
            X_series[pan]  = 1
            x_num += 1
        elif pan[0] == "X" and pan in X_series:
            X_series[pan] += 1
            x_num += 1
    
    X_all_series = {}
    for i in variants_all:
        try:
            pan = pango_lineage[i]
            if pan[0] == "X" and pan not in X_all_series:
                X_all_series[pan]  = 1
            elif pan[0] == "X" and pan in X_all_series:
                X_all_series[pan] += 1
            else:
                continue
        except:
            continue

    X_series_prop = {}
    for i in X_series:
        freq = round(X_series[i]/X_all_series[i],4)
        X_series_prop[i] = freq
    X_series_prop_sort = dict(sorted(X_series_prop.items(), key=lambda item:item[1],reverse=True))

    X_series_XY = {}
    X_series_XY_2 = {}
    for i in df.index:
        pan = df.loc[i,"pango_lineage"] 
        if pan[0] == "X" and pan not in X_series_XY:
            situ = df.loc[i,"lineage_X"] +", "+ df.loc[i,"lineage_Y"]
            X_series_XY[pan] = [situ]
        elif pan[0] == "X" and pan in X_series_XY:
            situ = df.loc[i,"lineage_X"] +", "+ df.loc[i,"lineage_Y"]
            X_series_XY[pan].append(situ)
            
    for i in X_series_XY:
        for j in set(X_series_XY[i]):
            if i not in X_series_XY_2:
                X_series_XY_2[i] = {}
                X_series_XY_2[i][j] = X_series_XY[i].count(j)
            else:
                 X_series_XY_2[i][j] = X_series_XY[i].count(j)
        X_series_XY_2[i] = dict(sorted(X_series_XY_2[i].items(), key=lambda item:item[1],reverse=True))
    
    with open("/home/soniali/Desktop/02_recom_230203/data/2023_02_08/3_recom_pattern/known_parental.txt") as f:
        rows =f.readlines()
    x_parents = {}
    for i in rows:
        x_parents[i.strip().split(":")[0]] = i.strip().split(":")[1].strip().split("[")[1].split("]")[0]

    table_result = {}
    for x in X_series:
        if x in x_parents:
            table_result[x] = X_all_series[x],X_series[x],X_series_prop_sort[x],x_parents[x],list(X_series_XY_2[x].keys())[0],list(X_series_XY_2[x].values())[0]
        else:
            table_result[x] = X_all_series[x],X_series[x],X_series_prop_sort[x],"  ",list(X_series_XY_2[x].keys())[0],list(X_series_XY_2[x].values())[0]

    df_identified_X_results = pd.DataFrame(columns = ["X_recombination_lineage","total recombinants","identified","identified rate", "known parents","identified parents","correspond number"],\
        index = range(len(table_result)))
    
    for i in df_identified_X_results.index:
        df_identified_X_results.loc[i,"X_recombination_lineage"] = list(table_result.keys())[i]
        for j in range(1,len(list(df_identified_X_results))):
            df_identified_X_results.loc[i,list(df_identified_X_results)[j]] = table_result[list(table_result.keys())[i]][j-1]
            
    df_identified_X_results = df_identified_X_results.sort_values(by = ["identified rate"],ascending = False)
    df_identified_X_results.to_csv("/home/soniali/Desktop/02_recom_230203/data/2023_02_08/3_recom_pattern/X_series_135567.csv",index = None)
    
    ## 
    df_identified_X_results = pd.read_csv("/home/soniali/Desktop/02_recom_230203/data/2023_02_08/3_recom_pattern/X_series_135567.csv")
    correspond_X_series = {}
    for i in range(46):
        if df_identified_X_results.loc[i,"comprehensive_consistency"] == "yes":
            correspond_X_series[df_identified_X_results.loc[i,"X_recombination_lineage"]] = df_identified_X_results.loc[i,"identified parents"]
    
    df_final = pd.read_csv(dirpath_pattern+"Independent_recombination_event_1451.csv")
    for i in df_final.index:
        epi = df_final.loc[i,"sample_id"]
        num_X = 0
        X_series = []
        if pango_lineage[epi][0] != "X":
            for e in each_group_sort[epi]:
                if pango_lineage[e][0] == "X":
                    num_X+=1
                    X_series.append(pango_lineage[e])
            if num_X >= 1 and (num_X/len(each_group_sort[epi])) >= 0.5:
                try:
                    recog = [df_final.loc[i,"lineage_X"],df_final.loc[i,"lineage_Y"]]
                    real = correspond_X_series[max(X_series, key = X_series.count)].split(", ")
                    if set(real) ==  set(recog):
                        df_final.loc[i,"X_series"] = max(X_series, key = X_series.count)
                    else:
                        df_final.loc[i,"X_series"] = "Others"
                except:
                    df_final.loc[i,"X_series"] = "Others"
            else:
                df_final.loc[i,"X_series"] = "Others"
        else:
            try:
                recog = [df_final.loc[i,"lineage_X"],df_final.loc[i,"lineage_Y"]]
                real = correspond_X_series[pango_lineage[epi]].split(", ")
                if set(real) ==  set(recog):
                    df_final.loc[i,"X_series"] = pango_lineage[epi]
                else:
                    df_final.loc[i,"X_series"] = "Others"
            except KeyError:
                df_final.loc[i,"X_series"] = "Others"

    # df_final.to_csv(dirpath_pattern+"Independent_recombination_event_"+str(df_final.shape[0])+"_v1.csv",index=None)
    ## Calculate the confidence value for each events
    # df_final = pd.read_csv(dirpath_pattern+"Independent_recombination_event_1451_v1.csv")
    all_lin = []
    for i in df_final.index:
        linX = del_star(df_final.loc[i,"lineage_X"])
        linY = del_star(df_final.loc[i,"lineage_Y"])
        all_lin.extend([linX,linY])

    for lin in set(all_lin):
        lin_name = lin.replace(".","")
        globals()["lin"+lin_name+"_epi"] = []

    for epi in pango_lineage:
        if set([pango_lineage[epi]]) - set(all_lin) == set():
            lin = pango_lineage[epi]
            lin_name = lin.replace(".","")
            globals()["lin"+lin_name+"_epi"].append(epi) 
        else:
            continue
    
    all_lin_name0 = []
    for i in df_final.index:
        linX, linY = df_final.loc[i,"lineage_X"], df_final.loc[i,"lineage_Y"]
        linX_name, linY_name = del_star(linX).replace(".",""), del_star(linY).replace(".","")
        all_lin_name0.extend([linX_name,linY_name])
    all_lin_name = set(all_lin_name0)
        
    lineage = {}
    for lin_name in all_lin_name:
        lin_epi = eval("lin"+lin_name+"_epi")
        record_mutY = {'count': 0, 'snps': {}}
        for seq_id in lin_epi:
            try:
                record_mutY['count'] += 1
                snps = variants_all[seq_id]
                for snp in snps:
                    record_mutY['snps'][snp] = record_mutY['snps'].get(snp, 0) + 1
            except:
                continue
        lineage[lin_name] = record_mutY
    
    for i in df_final.index:
        linX, linY = df_final.loc[i,"lineage_X"], df_final.loc[i,"lineage_Y"]
        X_muts, Y_muts = df_final.loc[i,"X_mutations"].split("/"), df_final.loc[i,"Y_mutations"].split("/")
        linX_name, linY_name = del_star(linX).replace(".",""), del_star(linY).replace(".","")
        Confidence_mean = get_linXY_mean_std(linX_name,lineage,Y_muts,linY_name,X_muts)

        df_final.loc[i,"Confidence_mean"] = Confidence_mean

    df_final.to_csv(dirpath_pattern+"Independent_recombination_event_"+str(df_final.shape[0])+"_v1.csv",index=None)
    data=df_final[['X_series', 'Confidence_mean']]
    data_X_series = data[data['X_series'] != 'Others']
    cultoff = round(np.percentile(np.array((sorted(data_X_series.Confidence_mean.tolist()))),25),4) #0.681
    df_final_2 = df_final.loc[df_final.Confidence_mean >= cultoff,]
    df_final_2.to_csv(dirpath_pattern+"Independent_recombination_event_1451_v1_hC939.csv",index=None)

if __name__ == '__main__':
    main()
    