__author__ = 'saltikov'

import pandas as pd
import heatmaps as hmp
import matplotlib.pyplot as plt


# The main file with RPKM values
df1 = pd.DataFrame.from_csv('/Users/saltikov/PycharmProjects/rna_seq/RUN1/RUN1_data.csv', index_col=0)

# columns with the specific RPKM values
cols = [0, 1, 2, 9, 10, 11, 18, 19, 20]

# Sorted list of DE genes
df2 = pd.DataFrame.from_csv('/Users/saltikov/PycharmProjects/rna_seq/RUN1/DESEQ/DESeq_RUN1_FCtable.csv',
                            index_col='Gene')

gs = hmp.GeneScan()
genes = gs.getgenes(df2, 'logFC_FumAs', 'p_FumAs', 2)
dfmap = gs.filterrpkm(df1, genes, cols)

print('{} genes survived.'.format(len(genes)))
print('The first 20 genes and their products:')
print(df1.ix[genes,['product']].head(20))
df1.ix[genes,['product']].to_clipboard()

hm = hmp.HeatMap(dfmap)
hm.heatmap_rpkm(dendro=True, labels=True)

coords = dict(title='/Users/saltikov/PycharmProjects/rna_seq/run1/RUN1_heatmap_FumVAs', output=True)
if coords['output']:
    plt.savefig(coords['title'] + '.png', format='png', dpi=200, bbox_inches='tight')

plt.show()