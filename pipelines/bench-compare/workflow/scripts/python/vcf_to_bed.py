import gzip


def lookup_maybe(a, b, k):
    return b[a.index(k)] if k in a else "."


with gzip.open(snakemake.input[0], "rt") as i, gzip.open(snakemake.output[0], "wt") as o:
    for line in i:
        if line.startswith("#"):
            pass
        elif ":FP:" in line or ":FN:" in line:
            s = line.rstrip().split("\t")
            chrom = s[0]
            start = int(s[1]) - 1
            ref = s[3]
            alt = s[4]
            end = start + len(ref)
            info = next((y[1] for x in s[7].split(";") if (y := x.split("="))[0] == "Regions"), ".")

            fmt = s[8].split(":")
            truth = s[9].split(":")
            query = s[10].split(":")

            reflen = len(ref)
            altlen = max([len(x) for x in alt.split(",")])

            if reflen == 1 and altlen == 1:
                vtype = "SNV"
            elif reflen >= 50 or altlen >= 50:
                vtype = "SV"
            elif reflen > 1 or altlen > 1 and reflen != altlen:
                vtype = "INDEL"
            else:
                vtype = "."

            newline = [
                chrom,
                str(start),
                str(end),
                ref,
                alt,
                info,
                lookup_maybe(fmt, truth, "BLT"),
                lookup_maybe(fmt, truth, "BD"),
                lookup_maybe(fmt, truth, "GT"),
                lookup_maybe(fmt, truth, "BK"),
                lookup_maybe(fmt, query, "BLT"),
                lookup_maybe(fmt, query, "BD"),
                lookup_maybe(fmt, query, "GT"),
                lookup_maybe(fmt, query, "BK"),
                str(reflen),
                str(altlen),
                vtype,
            ]

            o.write("\t".join(newline) + "\n")
