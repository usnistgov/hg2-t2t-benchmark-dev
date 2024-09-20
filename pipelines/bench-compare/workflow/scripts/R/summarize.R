library(tidyverse)

split_path <- function(df, path) {
  d <- basename(dirname(path))
  f <- basename(path)
  mutate(
    df,
    src_hap = str_extract(d, "hprc_([^_]+)_q100_([^_]+)", 1),
    dst_hap = str_extract(d, "hprc_([^_]+)_q100_([^_]+)", 2)
  )
}

gcols <- c("gchrom", "gstart", "gend", "realstart", "realend", "gid",
	   "overlap", "error_type", "q100", "hprc", "original",
	   "trim_left", "trim_right", "error_group", "error_group_size", 
	   "q100_len", "hprc_len", "seq"
	   )

gtypes = "ciiiiidcccciiiic"

vcols <- c("vchrom", "vstart", "vend", "vid", "ref", "alt", "regions", 
	   "truth_blt", "truth_bd", "truth_gt", "truth_bk",
	   "query_blt", "query_bd", "query_gt", "query_bk",
	   "truth_alt", "truth_len", 
	   "query_alt", "query_len", "genome_expected"
	   )

vtypes <- "ciiiccccccccccccicil"

read_vbench <- function(path) {
  read_tsv(
    path,
    col_types = c(vtypes, gtypes),
    col_names = c(vcols, gcols),
    na = ".",
  ) %>%
    split_path(path)
}

read_gbench <- function(path) {
  read_tsv(
    path,
    col_types = c(gtypes, vtypes),
    col_names = c(gcols, vcols),
    na = ".",
  ) %>%
    split_path(path)
}

snakemake@input[["gbench"]] %>%
  keep(~ str_detect(.x, "gbench")) %>%
  map_dfr(read_gbench) %>%
  write_tsv(snakemake@output[["gbench"]])

snakemake@input[["vbench"]] %>%
  keep(~ str_detect(.x, "vbench")) %>%
  map_dfr(read_vbench) %>%
  write_tsv(snakemake@output[["vbench"]])
