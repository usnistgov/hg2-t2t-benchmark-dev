from dataclasses import dataclass


# ASSUME sequences on the Q100 entry and HPRC entry are identical, so only keep
# one
@dataclass
class Error:
    #q100_chrom: str
    #q100_hap: str
    #q100_pos: int
    q100_seq: str
    hprc_contig: str
    hprc_pos: int
    hprc_seq: str
    hprc_dir: str


def split_error(q: str, h: str) -> Error:
    qs = q.split("_")
    hs = h.split("_")
    return Error(
        #q100_chrom=qs[0],
        #q100_hap=qs[1],
        #q100_pos=int(qs[2]),
        q100_seq=qs[3],
        hprc_seq=qs[4],
        hprc_contig=hs[0],
        hprc_pos=int(hs[1]),
        hprc_dir=hs[4],
    )


def star_maybe(s: str) -> str:
    return "*" if s == "" else s


def trim_error(e: Error) -> tuple[int, int]:
    # instances where Q100 is zero length seem to have a start position that is
    # off by one (to the left); fix that here
    if e.q100_seq == "*":
        return (1, 0)
    # do nothing for HPRC sequences that are zero length (nothing to trim)
    elif e.hprc_seq == "*":
        return (0, 0)
    # Otherwise trim both sequences as much as possible, starting from right
    else:
        ql = len(e.q100_seq)
        hl = len(e.hprc_seq)
        ml = min(ql, hl)
        i = 0
        while i < ml:
            if e.q100_seq[ql - i - 1] != e.hprc_seq[hl - i - 1]:
                break
            i = i + 1
        rm_right = i
        ml = ml - i
        i = 0
        while i < ml:
            if e.q100_seq[i] != e.hprc_seq[i]:
                break
            i = i + 1
        rm_left = i
        ql_right = ql - rm_right
        hl_right = hl - rm_right
        #e.q100_seq = star_maybe(e.q100_seq[rm_left:ql_right])
        #e.hprc_seq = star_maybe(e.hprc_seq[rm_left:hl_right])
        #e.hprc_pos = e.hprc_pos + rm_left
        return (0, 0)
        return (rm_left, rm_right)


# genome error ID
n = 0

with open(snakemake.input[0], "r") as ih, open(snakemake.output[0], "w") as oh:
    for i in ih:
        s = i.rstrip().split("\t")
        chrom = s[0]
        start = int(s[1])
        end = int(s[2])
        error_type = s[9]
        err = split_error(s[3], s[11])
        rm_l, rm_r = trim_error(err)
        newline = [
            chrom,
            str(start + rm_l),
            # add 1 to the end of zero-length variants since the
            # projection script needs regions >=1bp long
            str(end - rm_r + (1 if err.q100_seq == "*" else 0)),
            error_type,
            #err.q100_chrom,
            #err.q100_hap,
            #str(err.q100_pos),
            err.q100_seq,
            err.hprc_seq,
            f"{err.hprc_contig}_{err.hprc_pos}_{err.hprc_dir}",
            str(rm_l),
            str(rm_r),
            str(n),
        ]
        oh.write("\t".join(newline) + "\n")
        n = n + 1

