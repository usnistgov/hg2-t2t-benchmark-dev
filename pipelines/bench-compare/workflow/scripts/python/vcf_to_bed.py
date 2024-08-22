import gzip
import re


ip = snakemake.input[0]
op = snakemake.output[0]
hap = snakemake.wildcards["mapped_hap"]

if hap == "pat":
    is_pat = True
elif hap == "mat":
    is_pat = False
else:
    raise ValueError("womp womp")


def lookup_maybe(a, b, k):
    return b[a.index(k)] if k in a else "."


def select_alt(i, alts):
    if i == 0:
        return "."
    else:
        return alts[i-1]


def get_query_alt(gt, alt):
    g = re.match("([0-9])[/|]([0-9])", gt)
    if g is not None:
        a = alt.split(",")
        return (select_alt(int(g[1]), a), select_alt(int(g[2]), a))
    else:
        return (".", ".")


def get_variant_type(tlen, qlen):
    if tlen == "." or qlen == ".":
        return "."
    else:
        if tlen == 1 and qlen == 1:
            return "SNV"
        elif tlen >= 50 or qlen >= 50:
            return "SV"
        elif tlen > 1 or qlen > 1 and tlen != qlen:
            return "INDEL"
        else:
            return "UNK"


def get_variant_meta(tlen, query):
    qlen = "." if query == "." else len(query)
    vtype = get_variant_type(tlen, qlen)
    return vtype, str(qlen)


def star_maybe(s):
    return "*" if s == "" else s


def trim_ref(start, ref, truth, query1, query2):
    n = 0
    for r in ref:
        try:
            if truth != "." and truth[n] != r:
                break
            if query1 != "." and query1[n] != r:
                break
            if query2 != "." and query2[n] != r:
                break
        except IndexError:
            break
        n = n + 1
    if truth != ".":
        truth = truth[n:]
    if query1 != ".":
        query1 = query1[n:]
    if query2 != ".":
        query2 = query2[n:]
    return (start + n, ref[n:], truth, query1, query2)


with gzip.open(ip, "rt") as i, gzip.open(op, "wt") as o:
    id = 0
    for line in i:
        if line.startswith("#"):
            pass
        else:
        #elif ":FP:" in line or ":FN:" in line:
            id = id + 1

            s = line.rstrip().split("\t")

            # coordinates and other easy stuff
            chrom = s[0]
            start = int(s[1]) - 1
            ref = s[3]
            alt = s[4]
            info = next((y[1] for x in s[7].split(";") if (y := x.split("="))[0] == "Regions"), ".")

            fmt = s[8].split(":")
            truth = s[9].split(":")
            query = s[10].split(":")

            # next, split the ALT field so that we are only left with data relevant
            # for the desired haplotype
            truth_blt = lookup_maybe(fmt, truth, "BLT")
            truth_gt = lookup_maybe(fmt, truth, "GT")
            query_gt = lookup_maybe(fmt, query, "GT")

            alts = alt.split(",")
            query_alts = get_query_alt(query_gt, alt)

            if truth_blt == "het":
                # if truth is het, either keep or discard depending on which hap we want

                # ASSUME GT will always be defined for this case
                gt = re.match("([0-9])[/|]([0-9])", truth_gt)
                if gt is None:
                    raise ValueError("truth gt could not be found for het")
                pi = int(gt[1])
                mi = int(gt[2])

                # filtering paternal haplotype, paternal truth is not alt
                if is_pat and pi > 0:
                    truth_alt = alts[pi-1]
                # filtering maternal haplotype, maternal truth is not alt
                elif not is_pat and mi > 0:
                    truth_alt = alts[mi-1]
                # if neither of the above, then the desired hap must be ref in which
                # case we don't want it
                else:
                    # TODO keep everything
                    continue

            elif truth_blt == "homalt":
                # if truth is hom, always keep the line but split the alt field 

                # ASSUME this is defined for this case
                gt = re.match("([0-9])[/|]([0-9])", truth_gt)
                if gt is None:
                    raise ValueError("truth gt could not be found for homalt")
                truth_alt = alts[int(gt[1])-1]

            else:
                # otherwise this line is a FP with either nocall/halfcall for truth
                truth_alt = "."

            tlen = "." if truth_alt == "." else len(truth_alt)

            start_, ref_, truth_alt_, query1_alt, query2_alt = trim_ref(
                start,
                ref,
                truth_alt,
                query_alts[0],
                query_alts[1],
            )

            end = start_ + len(ref_)

            newline = [
                chrom,
                str(start_),
                str(end),

                str(id),

                star_maybe(ref_),
                alt,
                info,

                truth_blt,
                lookup_maybe(fmt, truth, "BD"),
                truth_gt,
                lookup_maybe(fmt, truth, "BK"),

                lookup_maybe(fmt, query, "BLT"),
                lookup_maybe(fmt, query, "BD"),
                query_gt,
                lookup_maybe(fmt, query, "BK"),

                star_maybe(truth_alt_),
                str(tlen),

                # 4 columns for each (potential) query field:
                # alt allele, variant type (INDEL, SNV), and alt length
                star_maybe(query1_alt),
                *get_variant_meta(tlen, query_alts[0]),

                star_maybe(query2_alt),
                *get_variant_meta(tlen, query_alts[1]),
            ]

            o.write("\t".join(newline) + "\n")
