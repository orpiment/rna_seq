# This code is based on:
# http://code.activestate.com/recipes/578834-hierarchical-clustering-heatmap-python/

__author__ = 'saltikov'
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
from collections import defaultdict
import os


class Chdir(object):
    """ Temporarily changes the working directory"""
    def __init__(self, newPath):
        self.savedPath = os.getcwd()
        self.newPath = newPath
    def __enter__(self):
        """Use a with Chdir as statement"""
        os.chdir(self.newPath)
    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

class HeatMap:
    def __init__(self, matrix):
        """ Generates a heat map using RPKM or other expression values

        :type matrix: pandas datafram indexed with gene ids
        :type self: Makes a heatmap
        """
        self.matrix = matrix

    def simpleheatmap(self):
        """Make a simple heatmap without labels and dendograms

        :return: get back the figure element
        """
        D = self.matrix
        D = D.replace(0, 1)
        fig = plt.figure(figsize=(6, 6.5))
        rect1 = [0.2, 0.1, 0.8, 0.8]
        ax1 = fig.add_axes(rect1)
        normal = mpl.colors.LogNorm(vmin=min(D.min()), vmax=max(D.max()))
        hm1 = ax1.matshow(D, cmap=plt.cm.jet, norm=normal, aspect='auto', origin='lower')
        ax1.set_yticks(np.arange(D.shape[0]))
        ax1.set_yticklabels([])
        ax1.set_xticks(np.arange(D.columns.shape[0]))
        ax1.set_xticklabels([])
        plt.colorbar(hm1, shrink=.4, pad=.1, aspect=10)
        return hm1

    def heatmap_rpkm(self, labels=False, dendro=False):
        """Make a heatmap with RPKM values.  Data should be in a pandas dataframe

        :param labels:
        :param dendro:
        :return:

        """
        ###**** Set up the figure ****###
        fig = plt.figure(num=1, figsize=(8, 10), dpi=72)
        fig.canvas.set_window_title('HeatMap_RPKM')
        D = self.matrix
        # data points of 0 cannot be plotted via log scale, convert 0 to 1
        D = D.replace(0, 1)

        ###**** make the side dendrogram ****###
        d1 = dist.pdist(D)
        D1 = dist.squareform(d1)
        rect1 = [0.72, 0.2, 0.1, 0.65]
        ax1 = fig.add_axes(rect1)
        Y1 = sch.linkage(D1, method='single', metric='braycurtis')
        Z1 = sch.dendrogram(Y1, orientation='left')
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ind1 = sch.fcluster(Y1, 0.5 * max(Y1[:, 2]), 'distance')
        idx1 = Z1['leaves']
        D = D.iloc[idx1, :]
        ind1 = ind1[idx1]

        ###***** make top dendrogram *****####
        d2 = dist.pdist(D.transpose())
        D2 = dist.squareform(d2)
        rect2 = [0.3, 0.86, 0.4, 0.1]
        ax2 = fig.add_axes(rect2)
        ax2.set_title('Heatmap analysis of RPKM for RUN1 RNA Seq Exp')
        Y2 = sch.linkage(D2, method='single', metric='braycurtis')
        Z2 = sch.dendrogram(Y2, orientation='top')
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ind2 = sch.fcluster(Y2, 0.7 * max(Y2[:, 2]), 'distance')
        idx2 = Z2['leaves']
        D = D.iloc[:, idx2]
        ind2 = ind2[idx2]

        #### make the heatmap using matshow ####
        rect3 = [0.3, 0.2, 0.4, 0.65]
        ax3 = fig.add_axes(rect3)
        normal = mpl.colors.LogNorm(vmin=min(D.min()), vmax=max(D.max()))
        # I tested using pcolor and pcolormesh would work but you will need to strip the index out
        # matshow seems to work the best
        hm = ax3.matshow(D, cmap=plt.cm.jet, norm=normal, aspect='auto')
        ax3.yaxis.set_ticks([])
        ax3.yaxis.set_ticklabels([])
        if labels:
            ax3.yaxis.set_ticks(np.arange(D.shape[0]))
            ax3.yaxis.set_ticklabels(D.index, minor=False, fontsize=7)
        ax3.xaxis.tick_bottom()
        ax3.xaxis.set_ticks(np.arange(D.columns.shape[0]))
        ax3.xaxis.set_ticklabels(D.columns, rotation=90, minor=False, fontsize=12)

        ###**** make colorbar scale ****###
        rect4 = [0.88, 0.3, 0.02, 0.4]
        ax4 = fig.add_axes(rect4)
        cb = plt.colorbar(hm, cax=ax4)
        ax4.set_title('RPKM', fontsize='small')

        ###**** make verticle colorbar for heatmap ****###
        rect5 = [0.705, 0.2, 0.01, 0.65]
        ax5 = fig.add_axes(rect5)
        dr = np.array(ind1, dtype=int)
        dr.shape = (len(ind1), 1)
        cmap = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
        ax5.matshow(dr, aspect='auto', origin='lower', cmap=cmap)
        ax5.set_xticks([])
        ax5.set_yticks([])

        return ax1, ax2, ax3, ax4, ax5

    def heatmap_fc(self, labels=False, dendro=False):
        """Make heatmap with log2fold change... function not ready yet

        :param labels:
        :param dendro:
        """
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


class GeneScan:
    def __init__(self):
        self._genes = []

    def filtergenes(self, df, log2fc_col, pval, logfc):
        """Main function for filtering a DESeq table for differentially expressed genes up and down regulated

        :param df: pandas dataframe of the DESeq data
        :param log2fc_col: particular column of the log2 foldchange
        :param pval: single value for p value, default is 0.05
        :param logfc: the criteria for filtering DE genes, e.g. 1 is 2 fold different

        """
        for key, item in df[(df[log2fc_col] >= logfc) | (df[log2fc_col] <= -logfc) & (df[pval] <= 0.05)].iterrows():
            if "Shewana3_R" in key:
                continue
            self._genes.append(key)

    def getgenes(self, df, log2fc_col, pval, log2fc=1):
        """Get a filtered list of genes from DESeq data

        :param df: pandas dataframe of the DESeq file
        :param log2fc_col: column in df that has the log2fc value
        :param pval: specifiy the pval, 0.05 is the default
        :param log2fc: specify the log2fc cutoff
        :return: list of genes
        \n
        =====Overview of the function=====\n
            Data should be generated from DESeq although another RNA seq type program might work.\n
        The data should be converted into a pandas dataframe indexed with the gene names.  In this\n
        case the column with the gene locus tags... e.g. Shewana3_1409.  The @log2fc_col string is \n
        the column from @df that contains the DESeq data.\n

        """
        self.filtergenes(df, log2fc_col, pval, log2fc)
        return self._genes

    def filterrpkm(self, df, genes, rpkm):
        """Get rpkm values for a list of genes from a panadas dataframe

        :param df: pandas dataframe with rpkm values
        :param genes: list of genes in df to search for
        :param rpkm: columns containing the rpkm values
        :return: a new dataframe of only the rpkm for the genes of interest
        """
        return df.ix[genes, rpkm]


class GeneScanNopval:
    def __init__(self):
        self._genes = []

    def filtergenes(self, df, log2fc_col, logfc):
        """Main function for filtering a DESeq table for differentially expressed genes

        :param df: pandas dataframe of the DESeq data
        :param log2fc_col: particular column of the log2 foldchange
        :param logfc: the criteria for filtering DE genes, e.g. 1 is 2 fold different

        """
        for key, item in df[(df[log2fc_col] >= logfc) | (df[log2fc_col] <= -logfc)].iterrows():
            if "Shewana3_R" in key:
                continue
            self._genes.append(key)

    def getgenes(self, df, log2fc_col, log2fc=1):
        """Get a filtered list of genes from DESeq data

        :param df: pandas dataframe of the DESeq file
        :param log2fc_col: column in df that has the log2fc value
        :param pval: specifiy the pval, 0.05 is the default
        :param log2fc: specify the log2fc cutoff
        :return: list of genes
        \n
        =====Overview of the function=====\n
            Data should be generated from DESeq although another RNA seq type program might work.\n
        The data should be converted into a pandas dataframe indexed with the gene names.  In this\n
        case the column with the gene locus tags... e.g. Shewana3_1409.  The @log2fc_col string is \n
        the column from @df that contains the DESeq data.\n

        """
        self.filtergenes(df, log2fc_col, log2fc)
        return self._genes

    def filterrpkm(self, df, genes, rpkm):
        """Get rpkm values for a list of genes from a panadas dataframe

        :param df: pandas dataframe with rpkm values
        :param genes: list of genes in df to search for
        :param rpkm: columns containing the rpkm values
        :return: a new dataframe of only the rpkm for the genes of interest
        """
        return df.ix[genes, rpkm]


if __name__ == '__main__':
    basedir = '/Users/saltikov/PycharmProjects/rna_seq/run2/'
    with Chdir(basedir) as cd1:
        # The main file with RPKM values
        df1 = pd.read_excel('RUN2_data.xlsx', 'sheet1', index_col=0)

    # columns with the specific RPKM values
    cols = [0, 1, 2, 9, 10, 11, 18, 19, 20]

    # Sorted list of DE genes from Maverixs using DESeq
    DEfiles = {'AsCr': ['AsIIIvsO2Cr_deseq_expression_filtered.tsv',
                        'O2AsIIIvsO2Cr_log2FoldChange', 'padj_O2AsIIIvsO2Cr'],
               'O2Cr': ['O2vsO2Cr_deseq_expression_filtered.tsv',
                        'OxygenvsO2Cr_log2FoldChange', 'padj_OxygenvsO2Cr'],
               'O2As': ['O2vsO2AsIII_deseq_expression_filtered.tsv',
                        'OxygenvsO2AsIII_log2FoldChange', 'padj_OxygenvsO2AsIII']}
    # Start to get list of DE genes
    gs = GeneScan()

    # Make a table containing the common DE gens
    table = defaultdict(int)

    # Placeholder for the genes
    genes = []

    # Filter criteria for log2 fold change in DESeq data
    # 1 would be a 2-fold change
    log2FC = 3

    with Chdir(basedir) as cd1:
        # Go through each DESeq file in the DEfiles dict and get the DE genes
        for condition, items in DEfiles.iteritems():
            df2 = pd.DataFrame.from_csv(items[0], sep='\t', index_col='id')
            found = (gs.getgenes(df2, items[1], items[2], log2FC))
            print('Found {} genes in {}'.format(len(found), condition))
            for item in found:
                table[item] += 1
        for gene in table.iterkeys():
            genes.append(gene)

    print('There are {} genes remaining'.format(len(genes)))

    # Get RPKM values for only the DE genes
    dfmap = gs.filterrpkm(df1, genes, cols)
    print(df1.ix[genes, ['product']].head(30))
    df1.ix[genes, ['product']]  #.to_clipboard()

    # call the HeatMap class
    hmap = HeatMap(dfmap)
    hmap.heatmap_rpkm(labels=True)

    coords = dict(title='RUN2_DEHeatMaps', output=False)
    if coords['output']:
        plt.savefig(coords['title'] + '.png', format='png', dpi=144, bbox_inches='tight')
    plt.show()