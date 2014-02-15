# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy.cluster.hierarchy as sch

# read in the excel file with all the data
df1 = pd.read_excel('/Users/saltikov/PycharmProjects/rna_seq/RUN2_data.xlsx', 'sheet1', index_col=0)
#df2 = pd.DataFrame.from_csv('O2vsO2Cr_deseq_expression_filtered.tsv', sep='\t', index_col='id')
df2 = pd.DataFrame.from_csv('/Users/saltikov/PycharmProjects/rna_seq/AsIIIvsO2Cr_deseq_expression_filtered.tsv', sep='\t', index_col='id')

# <codecell>

# these columns contain the RPKM values 0-2 = Arsenate
# 9-11 is Chromate and 18-20 is oxygen
cols = [0,1,2,9,10,11,18,19,20]

#just chromate and oxygen
#cols = [9,10,11,18,19,20]

# <codecell>

df2.columns.tolist()

# <codecell>

def get_gene_list(title, logchange=2, out=False):
    _coords = { 'genes' : [ ], 'save' : title , 'output' : out }
    for item in df2[df2.O2AsIIIvsO2Cr_log2FoldChange > logchange].index:
        if "Shewana3_R" in item:
            continue
        _coords['genes'].append(item) 
        # get the RPKM values of the genes listed in coords['genes']
    return df1.ix[_coords['genes'],cols], _coords

dfmap, coords = get_gene_list('Chromate_Expressed', 2.1, True)

dfmap.shape[0]
print len(coords['genes'])

# <codecell>

# need the values only.  This removes the gene names
D = dfmap.values
# can't have zero values for LogNorm scale
D = D + 1

# Compute and plot the side dendrogram.

# left, bottom, w, h
rectangle1 = (0,0.2,0.2,0.69)
fig = plt.figure(figsize=(7,9))
ax1 = fig.add_axes(rectangle1)
Y = sch.linkage(D, method='centroid')
Z1 = sch.dendrogram(Y, orientation='right')
ax1.set_xticks([])
ax1.set_yticks([])

# need to transpose the array so you can sort by RPKM
Dt = np.transpose(D)

# Compute and plot the top dendrogram.
ax2 = fig.add_axes([0.4,0.9,0.4,0.1])
Y = sch.linkage(Dt, method='single')
Z2 = sch.dendrogram(Y, p=50)
ax2.set_xticks([])
ax2.set_yticks([])

# Plot heatmap distance matrix.
axmatrix = fig.add_axes([0.4,0.2,0.4,0.69])
idx1 = Z1['leaves']
idx2 = Z2['leaves']
D = D[idx1 ,:]
D = D[: , idx2]

# set the heatmap color
color = plt.cm.jet
cmap = plt.get_cmap(color)

# normalize the colors based on min and max RPKMs
norm = mpl.colors.LogNorm(vmin=D.min(), vmax=D.max(), clip=True)

# this makes the heatmap
im = axmatrix.pcolormesh(D, cmap=cmap, norm=norm, clip_on=True)

# sort label elements idx1 and idx2 have a cerain order
# that must be preserved.  Use those lists to find the gene names
xlabels = list(dfmap.columns[i].__str__() for i in idx2)
ylabels = list(dfmap.index[i].__str__() for i in idx1)
axmatrix.set_xticklabels(xlabels, rotation=90, minor=False)
axmatrix.set_xticks(np.arange(xlabels.__len__())+0.5, minor=False)
axmatrix.set_yticklabels(ylabels, fontsize='small', minor=False)
axmatrix.set_yticks(np.arange(ylabels.__len__())+0.5, minor=False)
axmatrix.annotate('chromate\nresistance\noperon', xy=(7.5,11), xytext=(0.85,0.4), 
            textcoords='axes fraction', horizontalalignment="center",
            arrowprops=dict(facecolor='red', shrink=0.05))
axmatrix.axhspan(8, 11, 0, 9, edgecolor='k', facecolor="None", linewidth=2,
                 linestyle="dotted", label="chr operon")

axmatrix.annotate('nucleoside\nsalvage', xy=(7.5,50), xytext=(0.85,0.85), 
            textcoords='axes fraction', horizontalalignment="center",
            arrowprops=dict(facecolor='red', shrink=0.05))
axmatrix.axhspan(50, 52, 0, 9, edgecolor='k', facecolor="None", linewidth=2,
                 linestyle="dotted", label="chr operon")
plt.ylim(0, dfmap.shape[0])

# Plot colorbar.
axcolor = fig.add_axes([.85,0.2,0.02,0.6])
axcolor.set_title('RPKM')
plt.colorbar(im, cax=axcolor)


#plt.text(0.7, 1.1, coords['save'], horizontalalignment='center', transform = ax1.transAxes)

if coords['output']:
    plt.savefig(coords['save']+".png", format='png', dpi=200, bbox_inches='tight')

plt.show()

