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


def lookup_alt(fmt, sample, alts):
    """Return the alt for a given sample column and desired haplotype."""
    alt = "."
    gt = lookup_maybe(fmt, sample, "GT")
    # treat halfcalls as null
    if not (gt == "0|." or gt == ".|0"):
        # ASSUME anything unphased is a hom (in which case it doesn't matter)
        m = re.match("([0-9])[/|]([0-9])", gt)
        if m is not None:
            i = int(m[1 if is_pat else 2])
            if i > 0:
                alt = alts[i-1]
    return (gt, alt)


def star_maybe(s):
    return "*" if s == "" else s


def trim_alt_inner(s, n):
    if s == ".":
        return (".", ".")
    else:
        q = s[n:]
        return (star_maybe(q), str(len(q)))


def trim_ref(start, ref, truth, query):
    """Shave off extra bases on the front of ref and alts."""
    n = 0
    for r in ref:
        try:
            if truth != "." and truth[n] != r:
                break
            if query != "." and query[n] != r:
                break
        except IndexError:
            break
        n = n + 1
    _ref = ref[n:]
    newstart = start + n
    end = newstart + len(_ref)
    return (
        str(newstart), 
        str(end),
        star_maybe(_ref),
        trim_alt_inner(truth, n),
        trim_alt_inner(query, n),
    )


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

            # Split the alt field apart based on the haplotype we want.
            # ASSUME - ALL VARIANTS ARE PHASED; the index in the GT field
            # will be selected naively based on if it either the first
            # or second number (pat or mat resp.). A slash is allowed
            # but this is assumed to only apply to hom GT fields if present
            # at all.
            truth_blt = lookup_maybe(fmt, truth, "BLT")
            alts = alt.split(",")

            truth_gt, truth_alt = lookup_alt(fmt, truth, alts)
            query_gt, query_alt = lookup_alt(fmt, query, alts)
            
            # Skip if neither allele is an alt, which might happen if both are ref
            # and/or a halfcall or null value
            if truth_alt == "." and query_alt == ".":
                continue

            # Shave off prefix bases that are common to the ref/truth/query
            start_, end_, ref_, truth_, query_  = trim_ref(
                start,
                ref,
                truth_alt,
                query_alt,
            )

            newline = [
                chrom,
                start_,
                end_,

                str(id),

                ref_,
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

                *truth_,
                *query_,
                
            ]

            o.write("\t".join(newline) + "\n")
