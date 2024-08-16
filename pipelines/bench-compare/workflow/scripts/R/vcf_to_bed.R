library(tidyverse)

lookup_maybe <- function(x, k) {
  if (k %in% names(x)) { x[[k]] } else { NA }
}

split_format <- function(format, sample) {
  x <- str_split_1(sample, ":")
  names(x) <- str_split_1(format, ":")
  tibble(
    blt = lookup_maybe(x, "BLT"),
    bd = lookup_maybe(x, "BD"),
    gt = lookup_maybe(x, "GT"),
    bk = lookup_maybe(x, "BK"),
  )
}

snakemake@input[[1]] %>%
  read_tsv(
    col_types = "ci-cc---ccc",
    col_names = c("chrom", "start", "ref", "alt", "format", "truth", "query")
  ) %>%
  mutate(start = start - 1, end = str_length(ref) + start) %>%
  relocate(chrom, start, end) %>%
  mutate(truth = map2(format, truth, split_format)) %>%
  unnest(truth, names_sep = "_") %>%
  mutate(query = map2(format, query, split_format)) %>%
  unnest(query, names_sep = "_") %>%
  select(-format) %>%
  write_tsv(snakemake@output[[1]])
