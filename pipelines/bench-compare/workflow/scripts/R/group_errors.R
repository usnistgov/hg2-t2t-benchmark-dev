library(tidyverse)

chrom_index <- function(s) {
  s[s == "chrX"] <- "chr23"
  s[s == "chrY"] <- "chr24"
  as.integer(str_replace(s, "chr", ""))
}

nocov_ids <- read_tsv(
  snakemake@input[["nocov"]],
  col_names = "gid",
  col_types = "i"
) %>% pull(gid)

read_tsv(
  snakemake@input[["bed"]],
  col_types = "ciicccciii",
  col_names = c("chrom", "start", "end", "error_type", "q100", "hprc", "src", "trim_left", "trim_right", "gid")
) %>%
  separate(chrom, "_", into = c("chrom", "hap")) %>%
  mutate(chromidx = chrom_index(chrom)) %>%
  arrange(hap, chrom, start, end) %>%
  group_by(chrom, hap) %>%
  mutate(error_group = cumsum(start >= lag(cummax(end), default = TRUE))) %>%
  group_by(hap, chrom, error_group) %>%
  mutate(errorgroup_size = n()) %>%
  ungroup() %>%
  mutate(chrom = sprintf("%s_%s", chrom, hap)) %>%
  select(-chromidx, -hap) %>%
  mutate(nocov = gid %in% nocov_ids) %>%
  write_tsv(snakemake@output[[1]], col_names = F)
