library(tidyverse)

chrom_index <- function(s) {
  s[s == "chrX"] <- "chr23"
  s[s == "chrY"] <- "chr24"
  as.integer(str_replace(s, "chr", ""))
}

proc_df <- function(df) {
  df %>%
    separate_longer_delim(data, "~") %>%
    separate(data, ";", into = c(NA, NA, NA, NA, NA, NA, "gid", "error_group", "error_group_size", "nocov")) %>%
    mutate(gid = as.integer(gid)) %>%
    mutate(error_group_size = as.integer(error_group_size)) %>%
    mutate(chromidx = str_replace(chrom, "_[^_]+", "") %>% chrom_index()) %>%
    mutate(chrom = fct_reorder(factor(chrom), chromidx))
}

df0 <- read_tsv(
  snakemake@input[["src"]],
  col_names = c("chrom", "start", "end", "data"),
  col_types = "ciic"
) %>%
  proc_df

df1 <- read_tsv(
  snakemake@input[["projectable"]],
  col_names = c("chrom", "start", "end", "overlap", "data"),
  col_types = "ciidc"
) %>%
  proc_df

df2 <- read_tsv(
  snakemake@input[["smallvar"]],
  col_types = "c--iii----------",
  col_names = c("chrom", "start", "end", "gid")
)

small_ids <- df2$gid

df1_n <- df1 %>%
  group_by(gid) %>%
  tally()

df0 %>%
  left_join(df1_n, by = "gid") %>%
  replace_na(list(n = 0)) %>%
  ggplot(aes(x = n)) +
  geom_bar() +
  facet_wrap(c("chromidx")) +
  scale_y_log10() +
  labs(
    x = "copies on hg38",
    y = "count of unique errors"
  )
ggsave(snakemake@output[["counts_plot"]])

df0_grouped <- df0 %>%
  left_join(df1_n, by = "gid") %>%
  mutate(projected = !is.na(n)) %>%
  mutate(small = gid %in% small_ids) %>%
  mutate(
    group = case_when(
      projected & small ~ "small-variant",
      projected & !small ~ "projected",
      !projected & !small ~ "not-projected",
      TRUE ~ "unknown"
    )
  ) %>%
  replace_na(list(n = 0))

df0_grouped %>%
  ggplot(aes(x = factor(chromidx), fill = group)) +
  geom_bar() +
  labs(x = "chromosome", y = "unique error count")
ggsave(snakemake@output[["groups_plot"]])

df0_grouped %>%
  filter(small) %>%
  ggplot(aes(x = n)) +
  geom_bar() +
  facet_wrap(c("chromidx")) +
  scale_y_log10() +
  labs(
    x = "copies on hg38",
    y = "count of unique errors (smallvar only)"
  )
ggsave(snakemake@output[["smallvar_plot"]])

df0_grouped %>%
  filter(small) %>%
  group_by(chromidx, error_group_size) %>%
  tally() %>%
  ggplot(aes(error_group_size, n)) +
  geom_col() +
  facet_wrap(c("chromidx")) +
  scale_y_log10() +
  geom_label(aes(label = n), size = 3, position = position_dodge(width = 0.9)) +
  labs(x = "Overlapping set size", y = "unique error count (smallvar only)")
ggsave(snakemake@output[["smallvar_overlaps_plot"]])

df0_grouped %>%
  mutate(name = sprintf("%d_%d_%s", gid, n, group)) %>%
  mutate(strand = ".", score = 1000, thickStart = start, thickEnd = end) %>%
  mutate(itemRgb = case_when(small ~ "0,0,255", projected ~ "255,0,0", TRUE ~ "0,255,0")) %>%
  select(chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb, gid, n, group) %>%
  rename(hg38_copies = n) %>%
  write_tsv(snakemake@output[["groups_bed"]], col_names = F)

