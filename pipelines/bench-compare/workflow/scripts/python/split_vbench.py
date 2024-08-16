import gzip

ip = snakemake.input[0]
op = snakemake.output[0]

with gzip.open(ip, "rt") as i, gzip.open(op, "wt") as o:
    for line in i:
        pass
