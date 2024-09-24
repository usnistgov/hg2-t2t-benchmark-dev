library(tidyverse)

read_stuff <- function(path) {
  read_tsv(
    path,
    col_types = cols(
      regions = "c",
      nocov = "l",
      gstart = "i",
      gend = "i",
      genome_expected = "l",
      vid = "i",
      gid = "i",
      realstart = "i",
      realend = "i",
      overlap = "-",
      original = "-",
      q100_len = "i",
      hprc_len = "i",
      vstart = "i",
      vend = "i",
      truth_len = "i",
      query_len = "i",
      .default = "c",
      ),
    na = c(".", "NA")
  ) %>%
    mutate(
      liftover_diff = q100_len - (realend - realstart),
      gregion_len = realend - realstart,
      vregion_len = vend - vstart,
      ) %>%
    filter(src_hap == dst_hap) %>%
    select(-dst_hap) %>%
    rename(hap = src_hap) %>%
    mutate(
      across(
        c(ref, truth_alt, query_alt, q100, hprc),
        ~ if_else(.x == "*", "", .x)
      )
    ) %>%
    group_by(gid, hap) %>%
    mutate(ng = n()) %>%
    group_by(vid, hap) %>%
    mutate(nv = n()) %>%
    ungroup()
}

read_stuff(snakemake@input[[1]]) %>%
  filter(!gchrom %in% c("chrX", "chrY")) %>%
  mutate(gchrom = as.integer(str_sub(gchrom, 4))) %>%
  filter(abs(liftover_diff) < 100) %>%
  arrange(gid, desc(hap)) %>%
  filter(is.na(vchrom)) %>%
  select(-vchrom, -vstart, -vend, -vid, -ref, -alt, -regions, -starts_with("truth"),
         -starts_with("query"), -vregion_len, -nv, -ng, -genome_expected, -seq) %>%
  select(-gstart, -gend) %>%
  mutate(name = sprintf("%s_%s_%s_%s", hap, error_type, q100, hprc)) %>%
  mutate(gchrom = sprintf("chr%d", gchrom)) %>%
  relocate(gchrom, realstart, realend, name) %>%
  select(-error_type, -q100, -hprc, -q100_len, -hprc_len, -liftover_diff, -gregion_len) %>%
  write_tsv(snakemake@output[[1]], col_names = F)
