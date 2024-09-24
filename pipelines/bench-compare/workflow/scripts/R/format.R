library(tidyverse)

# dissimilarity = (indels + mismatches)/(projection block length) * 100
# TODO semicolon was a bad choice of separator since INFO uses it

chr2int <- function(s) {
  .s <- str_extract(s, "chr([0-9XY]+)_?.*", group = 1)
  if_else(.s == "X", 23, if_else(.s == "Y", 24, suppressWarnings(as.numeric(.s))))
}

int2chr <- function(i, suffix) {
  c <- if_else(i == 23, "X", if_else(i == 24, "Y", as.character(i)))
  sprintf("chr%s_%s", c, suffix)
}

## data_fmt <- c("REF", "ALT", "dissimilarity", "regions", "query_phase", "bench_phase") %>%
##   sprintf("%s=%%s", .) %>%
##   str_c(collapse = ":")

bed_path <- snakemake@input[["bed"]]
genome <- read_tsv(snakemake@input[["genome"]], col_names = c("chrom", "chr_length"), col_types = "ci")

# chrom_suffix <- if_else(str_detect(bed_path, "pat"), "PATERNAL", "MATERNAL")

readr::read_tsv(bed_path,
                col_types = "ciidc",
                col_names = c("chrom", "start", "end", "dissimilarity", "data")
                ) %>%
  mutate(chromidx = chr2int(chrom)) %>%
  separate_longer_delim(data, "~") %>%
  separate_wider_delim(data,
                       delim = ";",
                       names = c("error_type", "q100", "hprc", "src",
                                 "trim_left", "trim_right", "gid", "error_group", "error_group_size", "nocov")
                       ) %>%
  arrange(chromidx, chrom, start, end) %>%
  select(-chromidx) %>%
  left_join(genome, by = "chrom") %>%
  mutate(
    # variants with "*" are "between bases" and thus have zero length
    q100_len = if_else(q100 == "*", 0, str_length(q100)),
    hprc_len = if_else(hprc == "*", 0, str_length(hprc)),
    realstart = start,
    # fix the end for zero-length variants, had one base added to the end to
    # make the projection work
    realend = if_else(q100 == "*", pmax(start, end - 1), end),
    # add "slop"; some variants are near the beginning of a chromosome
    start = pmax(realstart - 50, 0),
    end = pmin(realend + 50, chr_length)
  ) %>%
  select(-chr_length) %>%
  # remove "SVs"
  filter(!(q100_len >= 50 | hprc_len >= 50)) %>%
  relocate(chrom, start, end, realstart, realend, gid) %>%
  readr::write_tsv(snakemake@output[[1]], col_names = FALSE)
