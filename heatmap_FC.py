__author__ = 'saltikov'
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def make_heatmap_fc(matrix):
    D =  matrix.values
    vmax, vmin = max(matrix.max()), min(matrix.min())
    colormap = my_colormap()
    normal = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    xlabels = list(matrix.columns)
    fig = plt.figure(figsize=(7, 9))
    axmatrix = fig.add_axes([0.4, 0.2, 0.4, 0.69])
    im = axmatrix.pcolormesh(D, cmap=colormap, norm=normal)
    axmatrix.set_xticklabels(xlabels, fontsize='small', rotation=90, minor=False)
    axmatrix.set_xticks(np.arange(xlabels.__len__()) + 0.5, minor=False)
    plt.ylim(0, D.shape[0])
    axcolor = fig.add_axes([.85, 0.2, 0.02, 0.6])
    plt.colorbar(im, cax=axcolor)
    axcolor.set_title('FC')
    # return fig

def my_colormap():
    cdict1 = [(0, (0, 0, 1)), # blue
              (0.2, (0, 1, 0)), # green
              (0.5, (0, 0, 0)), # black
              (0.6, (1, 1, 0)), # yellow
              (0.9, (1.0, 0.5, 0)), # orange
              (1, (1.0, 0, 0)), # red
     ]
    return mpl.colors.LinearSegmentedColormap.from_list('mycolors', cdict1)

if __name__ == '__main__':
    df1 = pd.read_excel('RUN1_Deseq_all_logFC.xlsx', sheetname='sheet1', index_col=0)
    cols = [ 3, 4, 5 ]

    df1map = df1.ix[:, cols]
    print df1map.dtypes

    D = df1map.values
    vmax, vmin = max(df1map.max()), min(df1map.min())
    print('Max: %f Min: %f' % (vmax, vmin))

    # do some filtering
    dfAs = df1[df1.ix[:,5] > 2]
    print dfAs.ix[:,(3,4,5,12)]

    make_heatmap_fc(dfAs.ix[:,(3,4,5)])
    plt.show()