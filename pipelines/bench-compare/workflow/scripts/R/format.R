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

path <- snakemake@input[[1]]

chrom_suffix <- if_else(str_detect(path, "pat"), "PATERNAL", "MATERNAL")

readr::read_tsv(path,
                col_types = "ciidc",
                col_names = c("chrom", "start", "end", "dissimilarity", "data")
                ) %>%
  mutate(chromidx = chr2int(chrom)) %>%
  separate_longer_delim(data, "~") %>%
  separate_wider_delim(data,
                       delim = ";",
                       names = c("name", "score", "strand", "thickStart",
                                 "thickEnd", "itemRgb", "consensus", "variant",
                                 "altName")
                       ) %>%
  ## mutate(
  ##   regions = str_replace(regions, "Regions=", ""),
  ##   query_phase = str_extract(SAMPLE, "(.*?):", group = 1),
  ##   bench_phase = str_extract(TRUTH, "(.*?):", group = 1),
  ##   data = sprintf(data_fmt, REF, ALT, dissimilarity, regions, query_phase, bench_phase)
  ## ) %>%
  ## filter(str_detect(regions, "notinsegdup")) %>%
  ## select(chrom, start, end, data) %>%
  arrange(chromidx, chrom, start, end) %>%
  select(-chromidx) %>%
  readr::write_tsv(snakemake@output[[1]], col_names = FALSE)
