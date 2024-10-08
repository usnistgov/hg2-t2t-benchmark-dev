library(tidyverse)

chrom_index <- function(s) {
  s[s == "chrX"] <- "chr23"
  s[s == "chrY"] <- "chr24"
  as.integer(str_replace(s, "chr", ""))
}

vcols <- c("vchrom", "vstart", "vend", "ref", "alt", "regions",
           "truth_blt", "truth_bd", "truth_gt", "truth_bk",
           "query_blt", "query_bd", "query_gt", "query_bk",
           "truth_alt", "truth_len",
           "query_alt", "query_len"
           )

vtypes <- "ciiccccccccccccici"

read_vbench <- function(path) {
  read_tsv(
    path,
    col_types = vtypes,
    col_names = vcols,
    na = ".",
  ) %>%
  mutate(chrom_idx = chrom_index(vchrom)) %>%
  arrange(chrom_idx, vstart, vend) %>%
  mutate(id = row_number()) %>%
  relocate(vchrom, vstart, vend, id) %>%
  select(-chrom_idx)
}

df <- read_vbench(snakemake@input[[1]])

# These should not have genome hits since ref and alt are exactly the same.
# Remove these first (the vast majority)
df_same <- df %>%
  filter(truth_alt == query_alt)

df_notsame <- df %>%
  anti_join(df_same, by = "id")

# Next find positions where a variant is split between two lines with the same
# ref, which often happens when truth is het and query is homalt (or reverse).
# In these cases, it is highly likely that (given a haplotype) the truth/query
# on the first line are X/NA and on the next line are NA/X. These are actually
# the same, and thus should not have a genome hit.
df_same2x <- df_notsame %>%
  group_by(vchrom, vstart, vend, ref) %>%
  filter(n() == 2) %>%
  replace_na(list(truth_alt = ".", query_alt = ".")) %>%
  arrange(vchrom, vstart, truth_alt, query_alt) %>%
  filter(truth_alt == lead(query_alt) | truth_alt == lag(query_alt)) %>%
  filter(n() == 2) %>%
  ungroup()

df_notsame2x <- df_notsame %>%
  anti_join(df_same2x, by = "id")

# Next remove "half TP" variants, which have a TP in truth/query and nocall in
# the other. Assume that these are generally pairs of variants that are on
# different positions (hence different lines with a nocall) in truth and query,
# but result in the same sequence and thus should not get a genome hit.
#
# NOTE this could also include variants that are wrongly marked TP due to
# vcfeval not taking phasing into account, but these should be rare.
df_same_split_tp <- df_notsame2x %>%
  filter(
  (truth_bd == "TP" & query_blt == "nocall")
  | (query_bd == "TP" & truth_blt == "nocall")
  )

# The remainder should be variants for which a genome hit is expected. Note that
# this should include variants marked TP for both truth/query with mismatched
# phase.
expected_ids <- df_notsame2x %>%
  anti_join(df_same_split_tp, by = "id") %>%
  pull(id)

df %>%
  mutate(genome_expected = id %in% expected_ids) %>%
  write_tsv(snakemake@output[[1]], col_names = FALSE)

