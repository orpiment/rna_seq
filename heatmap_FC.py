__author__ = 'saltikov'
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import defaultdict
import numpy as np
import pylab


def make_heatmap_fc(matrix, vmax, vmin, title):
    D = matrix.values
    colormap = my_colormap()
    normal = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    xlabels = list(matrix.columns)
    ylabels = list(matrix.index)
    fig = plt.figure(figsize=(4, 9))
    axmatrix = fig.add_axes([0.4, 0.2, 0.4, 0.69])
    im = axmatrix.pcolormesh(D, cmap=colormap, norm=normal)
    axmatrix.set_xticklabels(xlabels, fontsize='small', rotation=90, minor=False)
    axmatrix.set_xticks(np.arange(xlabels.__len__()) + 0.5, minor=False)
    axmatrix.set_title(title)
    if D.shape[0] < 80:
        axmatrix.set_yticklabels(ylabels, fontsize='x-small', minor=False)
        axmatrix.set_yticks(np.arange(ylabels.__len__()) + 0.5, minor=False)
    else:
        axmatrix.set_yticklabels([])
    plt.ylim(0, D.shape[0])
    axcolor = fig.add_axes([.85, 0.3, 0.03, 0.4])
    cbar = plt.colorbar(im, cax=axcolor)
    cbar.set_label('log Fold Change', fontsize='x-small')
    cbar.ax.tick_params(labelsize='x-small')
    axcolor.set_title('FC', fontsize='x-small')

def make_heatmap_fc2(matrix, vmax, vmin, title):
    D = matrix.values
    v = np.linspace(-8, 10.0, 15, endpoint=True)
    xlabels = list(matrix.columns)
    ylabels = list(matrix.index)
    fig = plt.figure(figsize=(4, 9))
    axmatrix = fig.add_axes([0.4, 0.2, 0.4, 0.69])
    im = axmatrix.pcolormesh(D, cmap=plt.cm.jet)
    axmatrix.set_xticklabels(xlabels, fontsize='small', rotation=90, minor=False)
    axmatrix.set_xticks(np.arange(xlabels.__len__()) + 0.5, minor=False)
    axmatrix.set_title(title)
    if D.shape[0] < 80:
        axmatrix.set_yticklabels(ylabels, fontsize='x-small', minor=False)
        axmatrix.set_yticks(np.arange(ylabels.__len__()) + 0.5, minor=False)
    else:
        axmatrix.set_yticklabels([])
    plt.ylim(0, D.shape[0])
    axcolor = fig.add_axes([.85, 0.3, 0.03, 0.4])
    cbar = plt.colorbar(im, ticks=v, cax=axcolor)
    cbar.set_label('log Fold Change', fontsize='x-small')
    cbar.ax.tick_params(labelsize='x-small')
    axcolor.set_title('FC', fontsize='x-small')


def my_colormap():
    cdict1 = [(0, (0, 0, 1)),  # blue
              (0.2, (0, 1, 0)),  # green
              (0.5, (0, 0, 0)),  # black
              (0.6, (1, 1, 0)),  # yellow
              (0.9, (1.0, 0.5, 0)),  # orange
              (1, (1.0, 0, 0)),  # red
    ]
    return mpl.colors.LinearSegmentedColormap.from_list('mycolors', cdict1)


def metabolic_cat(matrix):
    filtered = defaultdict(dict)
    for key, cat in matrix["Mainrole"].iteritems():
        filtered[cat]
    return list(x.__str__() for x in filtered.iterkeys())


class PlotMetabolic:
    def __init__(self, matrix):
        self.matrix = matrix

    def filter_metabolic_cat(self):
        self.mainrole = []
        self.number = []
        matrixF = self.matrix.groupby(self.matrix['Mainrole']).groups
        for key, item in matrixF.iteritems():
            self.mainrole.append(key)
            self.number.append(len(item))
        return self.mainrole, self.number

    def plot_metabolic_cat(self, yval, xval, title):
        self.mainrole = yval
        self.number = xval
        fig, ax1 = plt.subplots(figsize=(11, 7))
        plt.subplots_adjust(left=0.4, right=0.95)
        fig.canvas.set_window_title('Metabolic Catagories')
        pos = np.arange(len(self.mainrole)) + 0.5
        ax1.barh(pos, self.number, align='center', height=0.5, color='m')
        ax1.axis([0, max(self.number) + (max(self.number) * 0.1), 0, len(self.mainrole)])
        pylab.yticks(pos, self.mainrole, fontsize='x-small')
        ax1.set_title('Metabolic categories for: ' + title)


if __name__ == '__main__':
    # This is a processed file with DESeq results and gene categories
    # merged togther into one file
    df1 = pd.read_excel('RUN1_Deseq_all_logFC.xlsx', sheetname='sheet1', index_col=0)

    # these are the logFC values
    logFCcols = [ 3, 4, 5]

    # these are the p values
    pvalcols = [ 9, 10, 11]

    df1map = df1.ix[:, logFCcols]

    D = df1map.values
    vmax, vmin = max(df1map.max()), min(df1map.min())
    print('Max: %f Min: %f' % (vmax, vmin))

    metabolic_on = False
    heatmaps_on = True

    for x, y in zip(logFCcols, pvalcols):
        print x, y

        # do some filtering for logFC and p values
        dfRUN1 = df1[(df1.ix[:, x] < -1) | (df1.ix[:, x] > 1) & (df1.ix[:, y] < 0.05)]
        print(dfRUN1.ix[:, (3, 4, 5, 12)])

        if metabolic_on:
            mc = PlotMetabolic(dfRUN1)
            yval, xval = mc.filter_metabolic_cat()
            mc.plot_metabolic_cat(yval, xval, dfRUN1.columns[x])

        if heatmaps_on:
            make_heatmap_fc(dfRUN1.ix[:, logFCcols], vmax, vmin, dfRUN1.columns[x])
    plt.show()

