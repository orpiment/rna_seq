__author__ = 'saltikov'
import pandas as pd

basedir = '/Users/saltikov/Documents/Projects/RNA_Seq_Shewanella/RUN1/'

df1 = pd.read_csv(basedir+'RUN1_Deseq_all_logFC.txt', sep='\t', index_col=0)
dfGenesCoords = pd.DataFrame.from_csv(basedir+'gene_coordinates.txt', sep='\t', index_col='Gene')

infile = open(basedir+"endcoords.txt", 'r')
fileout = open(basedir+'locus_tags.txt', 'w')

for line in infile:
    line = int(line.rstrip())
    gene = dfGenesCoords[dfGenesCoords['stop_coord'] == line].index.tolist()
    print line, gene[0]
    fileout.write(gene[0]+'\n')
