__author__ = 'saltikov'
import heatmaps
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

df1 = pd.read_excel('RUN1_Deseq_all_logFC.xlsx', sheetname='sheet1', index_col=0)
cols = [ 3, 4, 5 ]

df1map = df1.ix[:, cols]

D = df1map.values

def make_heatmap1():
    fig = plt.figure(figsize=(7, 9))
    color = plt.cm.ocean
    colormap = plt.get_cmap(color)
    axmatrix = fig.add_axes([0.4, 0.2, 0.4, 0.69])
    normal = mpl.colors.Normalize(vmin=-10, vmax=10, clip=True)
    im = axmatrix.pcolormesh(D, cmap=colormap, norm=normal)
    plt.ylim(0, D.shape[0])
    axcolor = fig.add_axes([.85, 0.2, 0.02, 0.6])
    plt.colorbar(im, cax=axcolor)
    axcolor.set_title('FC')
    plt.show()

def make_heatmap2(matrix):
    D =  matrix.values
    fig = plt.figure(figsize=(7, 9))
    colormap = my_colormap()
    axmatrix = fig.add_axes([0.4, 0.2, 0.4, 0.69])
    normal = mpl.colors.Normalize(vmin=np.nanmin(D), vmax=np.nanmax(D))
    xlabels = list(df1.ix[:,cols].columns)
    im = axmatrix.pcolormesh(D, cmap=colormap, norm=normal)
    axmatrix.set_xticklabels(xlabels, fontsize='small', rotation=90, minor=False)
    axmatrix.set_xticks(np.arange(xlabels.__len__()) + 0.5, minor=False)
    plt.ylim(0, D.shape[0])
    axcolor = fig.add_axes([.85, 0.2, 0.02, 0.6])
    plt.colorbar(im, cax=axcolor).set_clim(-10,10)

    axcolor.set_title('FC')
    plt.show()

def my_colormap():
    cdict1 = [(0, (1, 1, 0)), # yellow
              (0.2, (1, 1, 0)), # yellow
              (0.3, (0, 1, 0)), # green
              (0.5, (0, 0, 0)), # black
              (0.65, (1.0, 0.65, 0)), # orange
              (0.73, (1.0, 0, 0)), # red
              (1, (1.0, 0, 0)), # red
     ]
    return mpl.colors.LinearSegmentedColormap.from_list('mycolors', cdict1)

dfAs = df1[df1.ix[:,5] < -4]
print dfAs.ix[:,(3,4,5,12)]
make_heatmap2(dfAs.ix[:,(3,4,5)])
#print(df1.head())