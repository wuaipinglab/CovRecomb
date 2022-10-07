

import matplotlib.pyplot as plt
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
    
    
dirpath = "/home/soniali/Desktop/03_CovRecomb/Simulation_Test/HOMO_probability/"
df_need_UA4 = pd.read_csv(dirpath+"HOMO_probability.txt")
df_draw = df_need_UA4[["mean_lin_diff","FPR","FDR","TPR","parallel_prop"]]

df_homo1 = df_draw[df_draw["parallel_prop"] == 0.1]
df_homo8 = df_draw[df_draw["parallel_prop"] == 0.5]
df_homo16 = df_draw[df_draw["parallel_prop"] == 0.8]

for i in [1,8,16]:
    df_temp = eval("df_homo"+str(i))
    df_temp = df_temp.sort_values(by = "mean_lin_diff",ascending=True)
    exec('fpr{} = {}'.format(i, list(df_temp["FPR"])))
    exec('tpr{} = {}'.format(i, list(df_temp["TPR"])))
    exec('fdr{} = {}'.format(i, list(df_temp["FDR"])))
    
    
name_list = [ int(i) for i in list(df_temp["mean_lin_diff"])]
x = list(range(len(name_list)))
total_width, n = 0.8, 3
width = total_width / n
bwith = 0.5  
line_color = "#c0c0c0"
plt.cla()
plt.figure(figsize=(8,3.2), dpi=300, facecolor='w')
plt.bar(x, fdr1, width=width, label='homo_pp0.1', tick_label=name_list, fc='#F0C8A3')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, fdr8, width=width, label='homo_pp0.5', fc='#7586AE')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, fdr16, width=width, label='homo_pp0.8', fc='#7C4D87')
for i in range(len(x)):
	x[i] = x[i] + width
plt.ylim((-0.01, 1.02))
plt.xlabel('Differential mutations between lineages', fontsize=8,fontweight = "bold")
plt.ylabel('False Discovery Rate', fontsize=8,fontweight = "bold")
font = {'weight': 'normal', 'size': 6}
plt.legend(prop=font)
y = MultipleLocator(0.2) 
ax = plt.gca()
ax.yaxis.set_major_locator(y)
bord_line(ax,bwith, line_color)
plt.savefig(dirpath+"FDR_UAB2_homopp.pdf")


plt.cla()
plt.figure(figsize=(8,3.2), dpi=300, facecolor='w')
plt.bar(x, tpr1, width=width, label='homo_pp0.1', tick_label=name_list, fc='#F0C8A3')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, tpr8, width=width, label='homo_pp0.5', fc='#7586AE')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, tpr16, width=width, label='homo_pp0.8', fc='#7C4D87')
for i in range(len(x)):
	x[i] = x[i] + width
plt.ylim((-0.01, 1.02))
plt.xlabel('Differential mutations between lineages', fontsize=8,fontweight = "bold")
plt.ylabel('True Discovery Rate', fontsize=8,fontweight = "bold")
font = {'weight': 'normal', 'size': 6}
plt.legend(prop=font)
y = MultipleLocator(0.2) 
ax = plt.gca()
ax.yaxis.set_major_locator(y)
bord_line(ax,bwith, line_color)
plt.savefig(dirpath+"TPR_UAB2_homopp.pdf")


plt.cla()
plt.figure(figsize=(8,3.2), dpi=300, facecolor='w')
plt.bar(x, fpr1, width=width, label='homo_pp0.1', tick_label=name_list, fc='#F0C8A3')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, fpr8, width=width, label='homo_pp0.5', fc='#7586AE')
for i in range(len(x)):
	x[i] = x[i] + width
plt.bar(x, fpr16, width=width, label='homo_pp0.8', fc='#7C4D87')
for i in range(len(x)):
	x[i] = x[i] + width
plt.ylim((-0.01, 1.02))
plt.xlabel('Differential mutations between lineages', fontsize=8,fontweight = "bold")
plt.ylabel('False Positive Rate', fontsize=8,fontweight = "bold")
font = {'weight': 'normal', 'size': 6}
plt.legend(prop=font)
y = MultipleLocator(0.2)
ax = plt.gca()
ax.yaxis.set_major_locator(y)
bord_line(ax,bwith, line_color)
plt.savefig(dirpath+"FPR_UAB2_homopp.pdf")
