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
  arrange(chromidx, chrom, start, end) %>%
  select(-chromidx, -score, -thickStart, -thickEnd, -itemRgb) %>%
  mutate(
    q100 = str_extract(name, "chr[0-9XY]+_[PM]ATERNAL_[0-9]+_([ACGT\\*]+)_[ACGT\\*]+", 1),
    hprc = str_extract(name, "chr[0-9XY]+_[PM]ATERNAL_[0-9]+_[ACGT\\*]+_([ACGT\\*]+)", 1),
    # variants with "*" are "between bases" and thus have zero length
    q100_len = if_else(q100 == "*", 0, str_length(q100)),
    hprc_len = if_else(hprc == "*", 0, str_length(hprc)),
    realstart = start,
    # for "*" variants on the q100 asm, shift end by 1, which usually will make it equal to
    # start unless the liftover did something weird
    realend = if_else(q100 == "*", pmax(end - 1, start), end),
    # some variants are near the beginning of a chromosome
    start = pmax(realend - 50, 0),
    end = realend + 50,
  ) %>%
  relocate(chrom, start, end, realstart, realend) %>%
  filter(!(q100_len >= 50 | hprc_len >= 50)) %>%
  readr::write_tsv(snakemake@output[[1]], col_names = FALSE)
