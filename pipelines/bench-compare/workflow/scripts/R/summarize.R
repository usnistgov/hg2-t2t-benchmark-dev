library(tidyverse)

split_path <- function(df, path, ext) {
  d <- basename(dirname(path))
  f <- basename(path)
  mutate(
    df,
    label = str_extract(f, "([^_]+)_.*", 1),
    membership = str_extract(f, sprintf(".*_([^_]+).%s", ext), 1),
    src_hap = str_extract(d, "hprc_([^_]+)_q100_([^_]+)", 1),
    dst_hap = str_extract(d, "hprc_([^_]+)_q100_([^_]+)", 2)
  )
}

read_vcf <- function(path) {
  read_tsv(
      path,
      col_types = "ci-cc---cc",
      col_names = c("chrom", "start", "ref", "alt", "format", "truth", "query")
    ) %>%
    split_path(path, "vcf.gz")
}

read_bed <- function(path) {
  read_tsv(
    path,
    "ciidc-----ccc",
    col_names = c("chrom", "start", "end", "overlap", "variant", "error_type", "var_type", "original")
  ) %>%
    split_path(path, "bed")
}

split_format <- function(format, sample) {
  x <- str_split_1(sample, ":")
  names(x) <- str_split_1(format, ":")
  tibble(blt = x[["BLT"]], bd = x[["BD"]], gt = x[["GT"]], bk = x[["BK"]])
}

snakemake@input[["gbench"]] %>%
  keep(~ str_detect(.x, "bed")) %>%
  map_dfr(read_bed) %>%
  arrange(chrom, start, end) %>%
  write_tsv(snakemake@output[["gbench"]])

snakemake@input[["vbench"]] %>%
  keep(~ str_detect(.x, "vcf")) %>%
  map_dfr(read_vcf) %>%
  mutate(
    reflen = str_length(ref),
    altlen = str_length(alt),
    var_type = case_when(
      reflen == 1 & altlen == 1 ~ "SNV",
      reflen != altlen & (reflen > 1 | altlen > 1) ~ "INDEL",
      # there should be none of these two
      reflen > 50 | altlen > 50 ~ "SV",
      TRUE ~ "UNK"
    )
  ) %>%
  arrange(chrom, start) %>%
  mutate(start = start - 1, end = start + reflen) %>%
  relocate(chrom, start, end) %>%
  mutate(truth = map2(format, truth, split_format)) %>%
  unnest(truth, names_sep = "_") %>%
  mutate(query = map2(format, query, split_format)) %>%
  unnest(query, names_sep = "_") %>%
  select(-format) %>%
  write_tsv(snakemake@output[["vbench"]])
