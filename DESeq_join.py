__author__ = 'saltikov'
import pandas as pd
import os

directory = '/Users/saltikov/Documents/Projects/RNA_Seq_Shewanella/RUN1/'
deseqdata = dict(FumAs=[], O2As=[], O2Fum=[])
extension = '.txt'
products = '/Users/saltikov/Documents/Projects/RNA_Seq_Shewanella/RUN1/Maverix_Files/products.txt'
categories = '/Users/saltikov/Documents/Projects/RNA_Seq_Shewanella/RUN1/Maverix_Files/Shewana3_genes_categories_v2.xlsx'

for key, item in deseqdata.iteritems():
    for f in os.listdir(directory+'DESEQ/'):
        if f.endswith(extension) and key in f:
            f = directory+'DESEQ/'+f
            item.append(f)

print deseqdata

def df_rename_combine(infile, name):
    df01 = pd.DataFrame()
    for i,file in enumerate(infile):
        df02 = pd.DataFrame.from_csv(file, sep='\s+', index_col='gene')
        df02.rename(columns = {'foldChange': 'FC_'+name+str(i),
                    'log2FoldChange' : 'log2FC_'+name+str(i),
                    'pval' : 'pval_'+name+str(i),
                    'padj': 'padj_'+name+str(i),
                    'baseMean' : 'baseMean_'+name+str(i),
                    'baseMeanA' : 'baseMeanA_'+name+str(i),
                    'baseMeanB' : 'baseMeanB_'+name+str(i),
                    },
                    inplace=True)
        df01 = df02.combine_first(df01)
    return df01

df1 = pd.DataFrame()
for key, item in deseqdata.iteritems():
    if key == 'FumAs':
        dfx0 = df_rename_combine(item, key)
        df1 = dfx0.combine_first(df1)
    if key == 'O2As':
        dfx1 = df_rename_combine(item, key)
        df1 = dfx1.combine_first(df1)
    if key == 'O2Fum':
        dfx2 = df_rename_combine(item, key)
        df1 = dfx2.combine_first(df1)

print('Getting products...')
df2 = pd.DataFrame.from_csv(products, sep='\t', index_col=0)

print('Joining products to main data frame...')
df3 = df1.join(df2)

print('Getting categories...')
df4 = pd.read_excel(categories, "sheet1", index_col=0)

print('Joining categories to main data frame...')
df5 = pd.merge(df3, df4, how='outer', left_index=True, right_index=True)

print('Making Excel file...')
df5.to_csv(directory+'RUN1_Deseq_all_logFC.txt', sep='\t')
#df5.to_excel(directory+'RUN1_Deseq_all_logFC.xlsx', sheet_name='sheet1')

print('Done!')