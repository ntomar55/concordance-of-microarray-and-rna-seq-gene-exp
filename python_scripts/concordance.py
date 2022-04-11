#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

drugs = ['FLUCONAZOLE', 'IFOSFAMIDE', 'LEFLUNOMIDE']
subsets = ['all', 'high', 'low']
iterables = [drugs, subsets]
idx = pd.MultiIndex.from_product(iterables, names=['Drug', 'subset'])
results = pd.DataFrame(index=idx, columns=['deseq', 'limma', 'concordance'])

for drug in drugs:
    deseq = pd.read_csv('/projectnb2/bf528/users/saxophone/data_p3/concordance/'
                        + drug + '_deseq_symbols.csv')
    limma = pd.read_csv('/projectnb2/bf528/users/saxophone/data_p3/concordance/'
                        + drug + '_limma_results_clean.csv')

    deseq.rename(columns = {'Unnamed: 0': 'symbol', 'log2FoldChange': 'logFC'},
                 inplace=True)
    limma.rename(columns = {'Unnamed: 0': 'probe', 'adj.P.Val': 'padj'},
                 inplace=True)

    print(drug)
    print(limma[:10][['symbol', 'probe', 'logFC', 'padj']])

    deseq_all = set(deseq.symbol)
    limma_all = set(limma.symbol)

    #print(drug)

    intersection = deseq_all & limma_all
    N = len(intersection)

    deseq = deseq[deseq.symbol.isin(intersection)]
    limma = limma[limma.symbol.isin(intersection)]

    def comp_concordance(deseq_set, limma_set, subset):
        deseq_diff = set(deseq_set[deseq_set.padj < 0.05].symbol)
        limma_diff = set(limma_set[limma_set.padj < 0.05].symbol)
        overlap = deseq_diff & limma_diff
        n0 = len(overlap)
        n1 = len(deseq_diff)
        n2 = len(limma_diff)
        nx = (N * n0 - n1 * n2) / (n0 + N - n1 - n2)
        concordance = 2 * nx / (n1 + n2)
        results.loc[(drug, subset)] = {
            'deseq': n1,
            'limma': n2,
            'concordance': concordance
        }

    comp_concordance(deseq, limma, 'all')

    median_exp = deseq.baseMean.median()
    deseq_high = deseq[deseq.baseMean > median_exp]
    limma_high = limma[limma.symbol.isin(deseq_high.symbol)]
    comp_concordance(deseq_high, limma_high, 'high')

    deseq_low = deseq[deseq.baseMean < median_exp]
    limma_low = limma[limma.symbol.isin(deseq_low.symbol)]
    comp_concordance(deseq_low, limma_low, 'low')

print(results)

to_graph = pd.DataFrame({subset: results.xs(subset, level='subset').concordance
                         for subset in subsets})
print(to_graph)

ax = to_graph.plot.bar(ylabel='Concordance')
ax.set_ylim(ymin=0)
plt.axes(ax)
#plt.axhline(0, color='gray')
plt.xticks(rotation='horizontal')
plt.savefig('/projectnb2/bf528/users/saxophone/data_p3/concordance_plot.png',
            dpi=300)
plt.show()

for subset in subsets:
    to_graph = results.xs(subset, level='subset')
    plt.scatter(to_graph.deseq, to_graph.concordance, label=subset)
    for drug in drugs:
        plt.annotate(drug[:3], (to_graph.loc[drug, 'deseq'],
                                to_graph.loc[drug, 'concordance']))

plt.legend(loc='lower right')
plt.xlabel('# D.E. genes from RNA-Seq')
plt.ylabel('Concordance')
plt.savefig('/projectnb2/bf528/users/saxophone/data_p3/seq_scatter.png',
            dpi=300)
plt.show()

for subset in subsets:
    to_graph = results.xs(subset, level='subset')
    plt.scatter(to_graph.limma, to_graph.concordance, label=subset)
    for drug in drugs:
        plt.annotate(drug[:3], (to_graph.loc[drug, 'limma'],
                                to_graph.loc[drug, 'concordance']))

plt.legend(loc='lower right')
plt.xlabel('# D.E. genes from microarray')
plt.ylabel('Concordance')
plt.savefig('/projectnb2/bf528/users/saxophone/data_p3/array_scatter.png',
            dpi=300)
plt.show()

