__author__ = 'saltikov'
import matplotlib.pylab as plt
import pandas as pd
import matplotlib as mpl
import numpy as np

def my_colormap():
    """Generates the color range for pcolormesh
    :return: matplotlib formated colorbar
    """
    cdict1 = [(0, (0, 0, 1)), # blue
              (0.2, (0, 1, 0)), # green
              (0.5, (0, 0, 0)), # black
              (0.6, (1, 1, 0)), # yellow
              (0.9, (1.0, 0.5, 0)), # orange
              (1, (1.0, 0, 0)), # red
    ]
    return mpl.colors.LinearSegmentedColormap.from_list('mycolors', cdict1)

# DESeq file containing logFoldChange and p values
df1 = pd.read_csv('DESeq_RUN1_FCtable.csv', index_col='Gene')
products = pd.read_csv('products.txt', sep='\t', index_col=0)

# columns in df1 that have the logFC values
cols = [ 1, 3, 5 ]

# make a new dataframe with logFC values
df1map = df1.ix[:, cols]

# select only highly expressed/repressed genes in a
# specific column of the dataframe
logFC = dict(coln=0, FC=[2.5, -2])
for item in logFC['FC']:
    if item > 0:
        dfinduced = df1map[df1map.iloc[:,logFC['coln']] > item]
    elif item < 0:
        dfrepressed = df1map[df1map.iloc[:,logFC['coln']] < item]

print('Joining categories to main data frame...')

# make a new dataframe with the hi values
dflogFC = dfinduced
dfup = pd.DataFrame('Up', index=dfinduced.index, columns=['As_v_Fum_Exp'])

dflogFC = dflogFC.append(dfrepressed)
dfdown = pd.DataFrame('Down', index=dfrepressed.index, columns=['As_v_Fum_Exp'])

# Make a list of the up/down genes with products for easy viewing
df2 = dfup.append(dfdown)
df3 = dflogFC.join(df2)
df4 = df3.join(products)
df4.to_csv('products_logFC.txt', sep='\t')

genes = dflogFC.index.tolist()

mymatrix, title = [ dfinduced, dfrepressed ], ['Induced', 'Repressed']
fig = plt.figure(figsize=(8, 7))
vmax, vmin = max(dflogFC.max()), min(dflogFC.min())
i = 0
for matrix in mymatrix:
    D = matrix.values
    colormap = my_colormap()
    normal = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    xlabels = list(matrix.columns)
    ylabels = list(matrix.index)
    print ylabels
    ax = plt.subplot(1,2,i+1)
    im = ax.pcolormesh(D, cmap=colormap, norm=normal, vmin=vmin, vmax=vmax)
    plt.title(title[i])
    plt.xticks(np.arange(xlabels.__len__()), xlabels, rotation=45, fontsize="x-small")
    plt.yticks(np.arange(ylabels.__len__())+0.5,ylabels, fontsize="x-small")
    plt.ylim(0, D.shape[0])
    cbar = plt.colorbar(im, shrink=0.6)
    cbar.set_label('log Fold Change', fontsize='x-small')
    i +=1

fig.subplots_adjust(left=0.13, right=0.96, bottom=0.15, wspace=0.7)
plt.savefig('heatmpas_logFC.png', dpi=200, bbox_inches='tight')
plt.show()

