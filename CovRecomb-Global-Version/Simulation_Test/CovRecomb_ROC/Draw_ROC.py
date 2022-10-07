
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import  auc
import pandas as pd
from matplotlib.ticker import MultipleLocator


def set_label(rects):
    for rect in rects:
        height = float(format(rect.get_height(),".2f"))
        plt.text(x = rect.get_x() + rect.get_width()/2,
                y = height + 0.05, 
                s = height, 
                ha = 'center',size = 5) 


def bord_line(ax,bwith, line_color):
    ax.spines['top'].set_color(line_color)
    ax.spines['right'].set_color(line_color)
    ax.spines['left'].set_color(line_color)
    ax.spines['bottom'].set_color(line_color)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['bottom'].set_linewidth(bwith)
    
    
dirpath = "/home/soniali/Desktop/03_CovRecomb/Simulation_Test/CovRecomb_ROC/"
df = pd.read_csv(dirpath+"ROC_results.txt")
df_need = df[df["seed_gen"] == 5]
df_need_UA1 = df_need[df_need["UAB"] == 1]
df_need_UA2 = df_need[df_need["UAB"] == 2]
df_need_UA3 = df_need[df_need["UAB"] == 3]
df_need_UA4 = df_need[df_need["UAB"] == 4]
# df_need_UA5 = df_need[df_need["UAB"] == 5]

for i in range(4):
    i += 1
    df_temp = eval("df_need_UA"+str(i))
    df_temp = df_temp.sort_values(by = "FPR",ascending=True)
    exec('fpr{} = {}'.format(i, list(df_temp["FPR"])))
    fpr_t = eval("fpr"+str(i))
    exec('tpr{} = {}'.format(i, list(df_temp["TPR"])))
    tpr_t = eval("tpr"+str(i))
    exec('auc{} = {}'.format(i, round(auc(fpr_t, tpr_t),3)))
    df_temp["mean_lin_diff"]
    exec('fdr{} = {}'.format(i, list(df_temp["FDR"])))

color_list = ["#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2","#BEB8DC","#E7DAD2","#999999"]
plt.cla()
plt.figure(figsize=(8, 7), dpi=300, facecolor='w')
plt.xlim((-0.01, 1.02))  
plt.ylim((-0.01, 1.02))
plt.xticks(np.arange(0, 1.1, 0.1))
plt.yticks(np.arange(0, 1.1, 0.1))
for i in range(4):
    i+=1
    plt.plot(eval("fpr"+str(i)), eval("tpr"+str(i)), color = color_list[i-1],lw=2, label=(i*"X")+(i*"Y"))
    
plt.legend(loc='lower right')
plt.xlabel('False Positive Rate', fontsize=14)  
plt.ylabel('True Positive Rate', fontsize=14)
plt.grid(b=True, ls=':')
plt.title(u'ROC curve', fontsize=18)
plt.savefig(dirpath+"ROC_UAB4.pdf")


df_temp = df_temp.sort_values(by = "mean_lin_diff",ascending=True)
name_list = [ int(i) for i in list(df_temp["mean_lin_diff"])]
x = list(range(len(name_list)))
total_width, n = 0.8, 4
width = total_width / n
bwith = 0.5  
line_color = "#c0c0c0"
plt.cla()
plt.figure(figsize=(8,3.2), dpi=300, facecolor='w')
plt.bar(x, fdr1, width=width, label='XY', tick_label=name_list, fc='#8ECFC9')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, fdr2, width=width, label='XXYY', fc='#FFBE7A')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, fdr3, width=width, label='XXXYYY', fc='#FA7F6F')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, fdr4, width=width, label='XXXXYYYY', fc='#82B0D2')
plt.ylim((-0.01, 1.02))
plt.xlabel('Differential mutations between lineages', fontsize=8,fontweight = "bold")
plt.ylabel('False Discovery Rate', fontsize=8,fontweight = "bold")
font = {'weight': 'normal', 'size': 6}
plt.legend(prop=font)
y = MultipleLocator(0.2) 
ax = plt.gca()
ax.yaxis.set_major_locator(y)
bord_line(ax,bwith, line_color)
plt.savefig(dirpath+"/FDR_UAB4.pdf")


df_temp = df_temp.sort_values(by = "mean_lin_diff",ascending=True)
name_list = [int(i) for i in list(df_temp["mean_lin_diff"])]
x = list(range(len(name_list)))
total_width, n = 0.8, 4
width = total_width / n
bwith = 0.5  
line_color = "#c0c0c0"
plt.cla()
plt.figure(figsize=(8,3.2), dpi=300, facecolor='w')
plt.bar(x, tpr1, width=width, label='XY', tick_label=name_list, fc='#8ECFC9')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, tpr2, width=width, label='XXYY', fc='#FFBE7A')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, tpr3, width=width, label='XXXYYY', fc='#FA7F6F')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, tpr4, width=width, label='XXXXYYYY', fc='#82B0D2')
plt.ylim((-0.01, 1.02))
plt.xlabel('Differential mutations between lineages', fontsize=8,fontweight = "bold")
plt.ylabel('True Discovery Rate', fontsize=8,fontweight = "bold")
font = {'weight': 'normal', 'size': 6}
plt.legend(prop=font)
y = MultipleLocator(0.2) 
ax = plt.gca()
ax.yaxis.set_major_locator(y)
bord_line(ax,bwith, line_color)
plt.savefig(dirpath+"/TPR_UAB4.pdf")


df_temp = df_temp.sort_values(by = "mean_lin_diff",ascending=True)
name_list = [ int(i) for i in list(df_temp["mean_lin_diff"])]
x = list(range(len(name_list)))
total_width, n = 0.8, 4
width = total_width / n
bwith = 0.5  
line_color = "#c0c0c0"
plt.cla()
plt.figure(figsize=(8,3.2), dpi=300, facecolor='w')
plt.bar(x, fpr1, width=width, label='XY', tick_label=name_list, fc='#8ECFC9')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, fpr2, width=width, label='XXYY', fc='#FFBE7A')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, fpr3, width=width, label='XXXYYY', fc='#FA7F6F')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, fpr4, width=width, label='XXXXYYYY', fc='#82B0D2')
plt.ylim((-0.01, 1.02))
plt.xlabel('Differential mutations between lineages', fontsize=8,fontweight = "bold")
plt.ylabel('False Positive Rate', fontsize=8,fontweight = "bold")
font = {'weight': 'normal', 'size': 6}
plt.legend(prop=font)
y = MultipleLocator(0.2) 
ax = plt.gca()
ax.yaxis.set_major_locator(y)
bord_line(ax,bwith, line_color)
plt.savefig(dirpath+"/FPR_UAB4.pdf")