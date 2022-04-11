#!/usr/bin/env python3

drugs = ['FLUCONAZOLE', 'IFOSFAMIDE', 'LEFLUNOMIDE']

for drug in drugs:
    deseq_all = set()
    deseq_diff = set()
    deseq_up = set()
    deseq_down = set()
    with open('/projectnb2/bf528/users/saxophone/data_p3/concordance/'
                   + drug + '_deseq_symbols.csv') as f:
        header = f.readline().rstrip('\n').replace('"', '').split(',')
        p_idx = -1
        l2fc_idx = -1
        for i in range(len(header)):
            if header[i] == 'padj':
                p_idx = i
            elif header[i] == 'log2FoldChange':
                l2fc_idx = i
        for line in f:
            fields = line.rstrip('\n').replace('"', '').split(',')
            gene = fields[0]
            p = fields[p_idx]
            deseq_all.add(gene)
            if p != 'NA' and float(p) < 0.05:
                deseq_diff.add(gene)
                if float(fields[l2fc_idx]) > 0:
                    deseq_up.add(gene)
                else:
                    deseq_down.add(gene)
    limma_all = set()
    limma_diff = set()
    limma_up = set()
    limma_down = set()
    with open('/projectnb2/bf528/users/saxophone/data_p3/concordance/'
                   + drug + '_limma_results_clean.csv') as f:
        header = f.readline().rstrip('\n').replace('"', '').split(',')
        p_idx = -1
        l2fc_idx = -1
        sym_idx = -1
        for i in range(len(header)):
            if header[i] == 'adj.P.Val':
                p_idx = i
            elif header[i] == 'logFC':
                l2fc_idx = i
            elif header[i] == 'symbol':
                sym_idx = i
        for line in f:
            fields = line.rstrip('\n').replace('"', '').split(',')
            gene = fields[sym_idx]
            p = float(fields[p_idx])
            limma_all.add(gene)
            if p < 0.05:
                limma_diff.add(gene)
                if float(fields[l2fc_idx]) > 0:
                    limma_up.add(gene)
                else:
                    limma_down.add(gene)

    print()
    print(drug)
    print('deseq results:')
    print('\t total # genes:', len(deseq_all))
    print('\t diff. expressed genes:', len(deseq_diff))
    print('\t upregulated:', len(deseq_up))
    print('\t downregulated:', len(deseq_down))
    print('limma results:')
    print('\t total # genes:', len(limma_all))
    print('\t diff. expressed genes:', len(limma_diff))
    print('\t upregulated:', len(limma_up))
    print('\t downregulated:', len(limma_down))

    print('Union of gene sets:', len(deseq_all | limma_all))
    print('Intersection of gene sets:', len(deseq_all & limma_all))
    print('In deseq, not in limma:', len(deseq_all - limma_all))
    print('In limma, not in deseq:', len(limma_all - deseq_all))
    print('Deseq significant, not in limma:', len(deseq_diff - limma_all))
    print('Limma significant, not in deseq:', len(limma_diff - deseq_all))
    continue
    input()
    for gene in deseq_diff - limma_all:
        print(gene)
    print('Limma sig, deseq absent:', len(limma_diff - deseq_all))
    input()
    for gene in  limma_diff - deseq_all:
        print(gene)
