#!/usr/bin/env python3

## Test which mappings in the provided file are unique.

from map_util import map_insert

#"REFSEQ","PROBEID","SYMBOL"

ref2probe = {}
ref2sym = {}
probe2ref = {}
probe2sym = {}
sym2ref = {}
sym2probe = {}


with open('/project/bf528/project_3/refseq_affy_map.csv') as f:
  f.readline()
  for line in f:
    line = line.replace('"', '').strip()
    fields = line.split(',')
    ref = fields[0]
    probe = fields[1]
    sym = fields[2]

    map_insert(ref, probe, ref2probe)
    map_insert(ref, sym, ref2sym)
    map_insert(probe, ref, probe2ref)
    map_insert(probe, sym, probe2sym)
    map_insert(sym, ref, sym2ref)
    map_insert(sym, probe, sym2probe)

for l in ref2probe.values():
  if len(l) > 1:
    print('REFSEQ-->PROBEID not unique')
    break

for l in ref2sym.values():
  if len(l) > 1:
    print('REFSEQ-->SYMBOL not unique')
    break

for k in probe2ref:
  if len(probe2ref[k]) > 1:
    print('PROBEID-->REFSEQ not unique')
    print(k, probe2ref[k])
    break

for k in probe2sym:
  if len(probe2sym[k]) > 1:
    print('PROBEID-->SYMBOL not unique')
    print(k, probe2sym[k])
    break

for l in sym2ref.values():
  if len(l) > 1:
    print('SYMBOL-->REFSEQ not unique')
    break

for l in sym2probe.values():
  if len(l) > 1:
    print('SYMBOL-->PROBEID not unique')
    break
