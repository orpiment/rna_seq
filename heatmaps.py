# <codecell>
__author__ = 'saltikov'
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy.cluster.hierarchy as sch

# <codecell>
class  HeatMap:
    def __init__(self, matrix):
        """ Generates a heat map using RPKM or other expression values

        :type matrix: pandas datafram indexed with gene ids
        :type self: Makes a heatmap
        """
        self.matrix = matrix

    def heatmap_rpkm(self, labels=False, dendro=False):
        fig = plt.figure(figsize=(7, 9))
        D = self.matrix.values
        D = D + 1
        if dendro:
            rectangle1 = (0, 0.2, 0.2, 0.69)
            ax1 = fig.add_axes(rectangle1)
            Y = sch.linkage(D, method='centroid')
            Z1 = sch.dendrogram(Y, orientation='right')
            ax1.set_xticks([])
            ax1.set_yticks([])

            # need to transpose the array so you can sort by RPKM
            Dt = np.transpose(D)

            # Compute and plot the top dendrogram.
            ax2 = fig.add_axes([0.4, 0.9, 0.4, 0.1])
            Y = sch.linkage(Dt, method='single')
            Z2 = sch.dendrogram(Y)
            ax2.set_xticks([])
            ax2.set_yticks([])

            # Plot heatmap distance matrix.
            axmatrix = fig.add_axes([0.4, 0.2, 0.4, 0.69])
            idx1 = Z1['leaves']
            idx2 = Z2['leaves']
            D = D[idx1, :]
            D = D[:, idx2]
        color = plt.cm.jet
        colormap = plt.get_cmap(color)
        axmatrix = fig.add_axes([0.4, 0.2, 0.4, 0.69])
        normal = mpl.colors.LogNorm(vmin=D.min(), vmax=D.max(), clip=True)
        im = axmatrix.pcolormesh(D, cmap=colormap, norm=normal, clip_on=True)
        if labels:
            if dendro:
                ylabels = list(self.matrix.index[i].__str__() for i in idx1)
            else:
                ylabels = list(self.matrix.index)
            axmatrix.set_yticklabels(ylabels, fontsize='small', minor=False)
            axmatrix.set_yticks(np.arange(ylabels.__len__()) + 0.5, minor=False)

        xlabels = list(self.matrix.columns[i].__str__() for i in idx2)
        axmatrix.set_xticklabels(xlabels, rotation=90, minor=False)
        axmatrix.set_xticks(np.arange(xlabels.__len__()) + 0.5, minor=False)

        plt.ylim(0, self.matrix.shape[0])
        axcolor = fig.add_axes([.85, 0.2, 0.02, 0.6])
        cbar = plt.colorbar(im, cax=axcolor)
        cbar.set_label('RPKM', fontsize='x-small')
        cbar.ax.tick_params(labelsize='x-small')

    def heatmap_fc(self, labels=False, dendro=False):
        fig = plt.figure(figsize=(7, 9))
        D = self.matrix.values
        if dendro:
            rectangle1 = (0, 0.2, 0.2, 0.69)
            ax1 = fig.add_axes(rectangle1)
            Y = sch.linkage(D, method='centroid')
            Z1 = sch.dendrogram(Y, orientation='right')
            ax1.set_xticks([])
            ax1.set_yticks([])

            # need to transpose the array so you can sort by RPKM
            Dt = np.transpose(D)

            # Compute and plot the top dendrogram.
            ax2 = fig.add_axes([0.4, 0.9, 0.4, 0.1])
            Y = sch.linkage(Dt, method='single')
            Z2 = sch.dendrogram(Y)
            ax2.set_xticks([])
            ax2.set_yticks([])

            # Plot heatmap distance matrix.
            axmatrix = fig.add_axes([0.4, 0.2, 0.4, 0.69])
            idx1 = Z1['leaves']
            idx2 = Z2['leaves']
            D = D[idx1, :]
            D = D[:, idx2]

        color = plt.cm.jet
        colormap = plt.get_cmap(color)
        axmatrix = fig.add_axes([0.4, 0.2, 0.4, 0.69])
        normal = mpl.colors.Normalize(vmin=np.nanmin(D), vmax=np.nanmax(D))
        im = axmatrix.pcolormesh(D, cmap=colormap, norm=normal, clip_on=True)
        if labels:
            if dendro:
                xlabels = list(self.matrix.columns[i].__str__() for i in idx2)
                ylabels = list(self.matrix.index[i].__str__() for i in idx1)
            else:
                ylabels = list(self.matrix.index)
                xlabels = list(self.matrix.columns)
            axmatrix.set_xticklabels(xlabels, rotation=90, minor=False)
            axmatrix.set_xticks(np.arange(xlabels.__len__()) + 0.5, minor=False)
            axmatrix.set_yticklabels(ylabels, fontsize='small', minor=False)
            axmatrix.set_yticks(np.arange(ylabels.__len__()) + 0.5, minor=False)
        plt.ylim(0, self.matrix.shape[0])
        axcolor = fig.add_axes([.85, 0.2, 0.02, 0.6])
        plt.colorbar(im, cax=axcolor)
        axcolor.set_title('FC')

class  GeneScan:
    def __init__(self):
        self._genes = []
    def filtergenes(self, df, log2fc, pval, fc=1.5):
        for key, item in df[(df[log2fc] >= fc) | (df[log2fc] <= -fc) & (df[pval] <= 0.05)].iterrows():
            if "Shewana3_R" in key:
                continue
            self._genes.append(key)
    def getgenes(self, df, log2fc, pval):
        self.filtergenes(df, log2fc, pval)
        return self._genes
    def filterrpkm(self, df, genes, rpkm):
        return df.ix[genes, rpkm]

def heatmap_dendro_rpkm(df, _coords):
    # need the values only.  This removes the gene names
    D = df.values

    # can't have zero values for LogNorm scale
    D = D + 1

    # make a figure instance
    fig = plt.figure(figsize=(7, 9))

    # Compute and plot the side dendrogram.
    # left, bottom, w, h
    rectangle1 = (0, 0.2, 0.2, 0.69)
    ax1 = fig.add_axes(rectangle1)
    Y = sch.linkage(D, method='centroid')
    Z1 = sch.dendrogram(Y, orientation='right')
    ax1.set_xticks([])
    ax1.set_yticks([])

    # need to transpose the array so you can sort by RPKM
    Dt = np.transpose(D)

    # Compute and plot the top dendrogram.
    ax2 = fig.add_axes([0.4, 0.9, 0.4, 0.1])
    Y = sch.linkage(Dt, method='single')
    Z2 = sch.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])

    # Plot heatmap distance matrix.
    axmatrix = fig.add_axes([0.4, 0.2, 0.4, 0.69])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = D[idx1, :]
    D = D[:, idx2]

    color = plt.cm.jet
    cmap = plt.get_cmap(color)
    #levels = mpl.ticker.LogLocator(base=10, subs=[1], numticks=20).tick_values(1, D.max())
    #norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    norm = mpl.colors.LogNorm(vmin=D.min(), vmax=D.max(), clip=True)

    im = axmatrix.pcolormesh(D, cmap=cmap, norm=norm, clip_on=True)

    # sort label elements
    xlabels = list(df.columns[i].__str__() for i in idx2)
    ylabels = list(df.index[i].__str__() for i in idx1)
    axmatrix.set_xticklabels(xlabels, rotation=90, minor=False)
    axmatrix.set_xticks(np.arange(xlabels.__len__()) + 0.5, minor=False)
    axmatrix.set_yticklabels(ylabels, fontsize='small', minor=False)
    axmatrix.set_yticks(np.arange(ylabels.__len__()) + 0.5, minor=False)
    plt.ylim(0, df.shape[0])

    # Plot colorbar.
    axcolor = fig.add_axes([.85, 0.2, 0.02, 0.6])
    plt.colorbar(im, cax=axcolor)
    axcolor.set_title('RPKM')

    if _coords['output']:
        plt.savefig(_coords['title'] + '.png', format='png', dpi=200, bbox_inches='tight')
    plt.show()

# <codecell>
if __name__ == '__main__':

    # The main file with RPKM values
    df1 = pd.read_excel('/Users/saltikov/PycharmProjects/rna_seq/RUN2_data.xlsx', 'sheet1', index_col=0)

    # columns with the specific RPKM values
    cols = [0, 1, 2, 9, 10, 11, 18, 19, 20]

    # Sorted list of DE genes
    df2 = pd.DataFrame.from_csv('AsIIIvsO2Cr_deseq_expression_filtered.tsv', sep='\t', index_col='id')

    gs = GeneScan()
    genes = gs.getgenes(df2, 'O2AsIIIvsO2Cr_log2FoldChange', 'padj_O2AsIIIvsO2Cr')
    dfmap = gs.filterrpkm(df1, genes, cols)

    print('{} genes survived.'.format(len(genes)))
    print('The first 20 genes and their products:')
    print(df1.ix[genes,['product']].head(20))

    # heatmap_dendro(dfmap, coords)

    hm = HeatMap(dfmap)
    hm.heatmap_rpkm(dendro=True, labels=False)
    plt.show()