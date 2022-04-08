#!/usr/bin/env python3

import map_util as mu

# Build map of refs to symbols
ref2symbol = mu.create_map('REFSEQ', 'SYMBOL')

conditions = ['CAR_PXR', 'DAN_Dmg', 'AhR']
drugs = ['FLUCONAZOLE', 'IFOSFAMIDE', 'LEFLUNOMIDE']

for i in range(3):
    cond = conditions[i]
    drug = drugs[i]
    # Avoid duplicates; use the most differentially expressed ref (will be
    # the first in the file, since file is sorted by p-value)
    num_amb = 0
    num_missing = 0
    num_duplicate = 0
    used_symbols = set()
    with open(('/projectnb2/bf528/users/saxophone/DESeq_results/'
               + cond + '_vs_ctrl.csv'), 'r') as infile:
        with open(('/projectnb2/bf528/users/saxophone/data_p3/concordance/'
                   + drug + '_deseq_symbols.csv'), 'w') as outfile:
            header = infile.readline()
            outfile.write(header)
            for line in infile:
                first_comma = line.find(',')
                ref = line[:first_comma].replace('"', '')
                if ref not in ref2symbol:
                    num_missing += 1
                    print("Can't find ref", ref)
                    input('press enter to continue')
                    continue
                symbol_list = ref2symbol[ref]
                if len(symbol_list) == 1:
                    symbol = symbol_list[0]
                else:
                    num_amb += 1
                    continue
                    print(f'Possible symbols for {ref}:')
                    for i in len(symbol_list):
                        print(i, symbol_list[i])
                    selection = -1
                    prompt = 'Select the number of the correct gene symbol: '
                    while selection < 0 or selection >= len(symbol_list):
                        select_text = input(prompt)
                        try:
                            selection = int(select_text)
                        except ValueError:
                            selection = -1
                if symbol not in used_symbols:
                    used_symbols.add(symbol)
                    outfile.write(symbol)
                    outfile.write(line[first_comma:])
                else:
                    num_duplicate += 1
    print('Stats for', drug)
    print('Number of ambiguous refs:', num_amb)
    print('Number of missing refs:', num_missing)
    print('Number of refs mapped to genes:', len(used_symbols))
    print('Number of duplicate refs:', num_duplicate)
