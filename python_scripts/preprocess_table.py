#!/usr/bin/env python3

## Preprocess supp_table_14.csv to remove first line and repeat gene name
## on each line.

from os import chdir

chdir('/projectnb2/bf528/users/saxophone/data_p3/')
probe2sym = {}
ref2sym = {}

with open('supp-table-14.tsv', encoding="utf-8") as in_file:
    with open('probe-mapping.csv', 'w') as probe_out, open('ref-mapping.csv', 'w') as ref_out:
        in_file.readline()
        header = in_file.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            if header[i] == 'RefSeq':
                symidx = i
            if header[i] == 'Affymetrix':
                prbidx = i
            if header[i] == 'P2':
                p2idx = i
            if header[i] == 'P3':
                p3idx = i
        linenum = 2
        for line in in_file:
            linenum += 1
            fields = line.rstrip('\n').replace('"', '').split('\t')
            if len(fields[symidx]) > 0:
                sym = fields[symidx]
                # if line is blank, corresponds to previous gene symbol
            probe = fields[prbidx]
            p2 = fields[p2idx]
            p3 = fields[p3idx]

            if probe in probe2sym:
                old = probe2sym[probe]
                if sym != old:
                    if sym in old:
                        probe2sym[probe] = sym
                    elif old not in sym:
                        print(f'Line {linenum}: Probe '
                                f'{probe} could map to either '
                                f'{old} or {sym}')
                        response = 'n'
                        while response != 'y':
                            in_sym = input('Enter correct mapping: ')
                            response = input(f'Is "{in_sym}" correct? y/n ')
                        sym = in_sym
            else:
                probe2sym[probe] = sym

            #for ref in [*p2.split(','), *p3.split(',')]:
            for ref in p3.split(','):
                if ref in ref2sym:
                    old = ref2sym[ref]
                    if sym != old:
                        if sym in old:
                            ref2sym[ref] = sym
                        elif old not in sym:
                            print(f'Line {linenum}: Reference sequence '
                                  f'{ref} could map to either '
                                  f'{old} or {sym}')
                            response = 'n'
                            while response != 'y':
                                in_sym = input('Enter correct mapping: ')
                                response = input(f'Is "{in_sym}" correct? y/n ')
                            sym = in_sym
                elif len(ref) > 0:
                    ref2sym[ref] = sym
for probe in probe2sym:
    if probe2sym[probe].find(',') >= 0:
        input(probe + ': ' + probe2sym[probe])
for ref in ref2sym:
    if ref2sym[ref].find(',') >= 0:
        input(ref + ': ' + ref2sym[ref])
