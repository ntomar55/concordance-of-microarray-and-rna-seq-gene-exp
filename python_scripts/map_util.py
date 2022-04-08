### Defines helper functions for mapping probes, reference sequences, and
### gene symbols.

def map_insert(k, v, d):
    if k in d:
        if v in d[k]:
            return
        else:
            d[k].append(v)
    else:
        d[k] = [v]

def create_map(key, value):
    """
    Reads `refseq_aff_map.csv` and constructs a dictionary from the column
    whose label matches `key` to the column whose label matchs `value`.
    Valid labels are: "REFSEQ", "PROBEID", "SYMBOL", "GENENAME"
    """
    with open('/project/bf528/project_3/refseq_affy_map.csv') as f:
        mapping = {}
        header = f.readline().replace('"', '').strip().split(',')
        for i in range(len(header)):
            if header[i] == key:
                kidx = i
            if header[i] == value:
                vidx = i
        for line in f:
            fields = line.replace('"', '').strip().split(',')
            k = fields[kidx]
            v = fields[vidx]
            map_insert(k, v, mapping)
        return mapping
