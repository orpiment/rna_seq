__author__ = 'saltikov'
import matplotlib.pyplot as plt
import pandas as pd
fig, (ax1, ax2) = plt.subplots(2,1)
fig.set_size_inches(10,12)
# open separate data files
input_files = "output.csv"
growth = pd.DataFrame.from_csv(input_files, index_col=0)
ax1.errorbar(growth.index,growth[['A1', 'A2', 'A3']].mean(axis=1), label='wt 5mM')
ax1.errorbar(growth.index,growth[['B1', 'B2', 'B3']].mean(axis=1), label='0058 5mM')
ax1.errorbar(growth.index,growth[['C1', 'C2', 'C3']].mean(axis=1), label='0537 5mM')
ax1.errorbar(growth.index,growth[['D1', 'D2', 'D3']].mean(axis=1), label='0857 5mM')

ax2.errorbar(growth.index,growth[['A1', 'A2', 'A3']].mean(axis=1), label='wt 5mM')
ax2.errorbar(growth.index,growth[['E1', 'E2', 'E3']].mean(axis=1), label='2062 5mM')
ax2.errorbar(growth.index,growth[['F1', 'F2', 'F3']].mean(axis=1), label='2512 5mM')
ax2.errorbar(growth.index,growth[['G1', 'G2', 'G3']].mean(axis=1), label='3226 5mM')
ax2.errorbar(growth.index,growth[['H1', 'H2', 'H3']].mean(axis=1), label='3478 5mM')

ax1.set_ylabel('Growth OD' r'$_{600nm}$')
ax1.set_xlabel('Time (min)')
ax2.set_ylabel('Growth OD' r'$_{600nm}$')
ax2.set_xlabel('Time (min)')

lines1, labels1 = ax1.get_legend_handles_labels()
ax1.legend(lines1, labels1, bbox_to_anchor=(1.3, 1), ncol=1, fancybox=True, shadow=True)

lines1, labels1 = ax2.get_legend_handles_labels()
ax2.legend(lines1, labels1, bbox_to_anchor=(1.3, 1), ncol=1, fancybox=True, shadow=True)


plt.show()