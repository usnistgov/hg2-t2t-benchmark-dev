library(tidyverse)

pretty_theme <-
  theme(text = element_text(size = 6),
        line = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.25, "mm"),
        axis.line = element_line(linewidth = 0.15),
        legend.box.spacing = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0.75, 0, "mm"),
        plot.margin = margin(1, 1, 1, 1, "mm"),
        legend.key.size = unit(0.55, "lines"),
        legend.spacing = unit(1.5, "mm"),
        strip.text = element_text(size = rel(1), margin = margin(1, 1, 1, 1, "mm")),
        strip.background = element_rect(linetype = "blank")
        )
        

df_happy <- read_csv(
  snakemake@input[["happy"]],
  col_types = cols(
    Type = "c",
    Subtype = "c",
    Subset = "c",
    Filter = "c",
    Genotype = "-",
    QQ.Field = "-",
    QQ = "-",
    TRUTH.TOTAL.ti = "c",
    TRUTH.TOTAL.tv = "c",
    .default = "d"
  ),
  na = c(".", "")
) %>%
  filter(Filter == "PASS") %>%
  select(-Filter) %>%
  select(-matches("\\.(ti|tv)"))

important_subsets <- list(
  c("XY", "AllAutosomes"),
  c("LowComplexity", "AllTandemRepeatsandHomopolymers_slop5"),
  ## c("LowComplexity", "AllHomopolymers_ge7bp_imperfectge11bp_slop5"),
  c("LowComplexity", "SimpleRepeat_diTR_10to49_slop5"),
  c("LowComplexity", "SimpleRepeat_homopolymer_ge21_slop5"),
  c("LowComplexity", "SimpleRepeat_homopolymer_ge12_slop5"),
  c("LowComplexity", "SimpleRepeat_imperfecthomopolge11_slop5"),
  c("LowComplexity", "SimpleRepeat_imperfecthomopolge21_slop5"),
  c("LowComplexity", "SimpleRepeat_quadTR_50to149_slop5"),
  c("LowComplexity", "SimpleRepeat_quadTR_50to149_slop5"),
  c("LowComplexity", "SimpleRepeat_triTR_50to149_slop5"),
  c("Mappability", "lowmappabilityall"),
  c("Union", "alldifficultregions"),
  ## c("Union", "alllowmapandsegdupregions"),
  c("GCcontent", "gclt25orgt65_slop50"),
  c("SegmentalDuplications", "segdups")
)

phred <- function(x) {
  - 10 * log10(x)
}

read_error_bed <- function(path, ispat) {
  read_tsv(
    path,
    col_types = "c--ii--c--------ii-",
    col_names = c("chrom", "start", "end", "error_type", "q100_len", "hprc_len")
  ) %>%
    add_column(ispat = ispat) %>%
    filter(!chrom %in% c("chrX", "chrY")) %>%
    mutate(chrom = as.integer(str_sub(chrom, 4))) %>%
    mutate(Type = if_else(q100_len == 1 & hprc_len == 1, "SNP", "INDEL"))
}

df_error_bed <- bind_rows(
  read_error_bed(snakemake@input[["mat_errors"]], FALSE),
  read_error_bed(snakemake@input[["pat_errors"]], TRUE)
)

strat_root <- snakemake@input[["all_strats"]]

read_strat <- function(level, name) {
  read_tsv(
    file.path(strat_root, level, sprintf("GRCh38_%s.bed.gz", name)),
    col_types = "cii",
    col_names = c("s_chrom", "s_start", "s_end")
  ) %>%
    add_column(level = level, name = name) %>%
    filter(!s_chrom %in% c("chrX", "chrY")) %>%
    mutate(s_chrom = as.integer(str_sub(s_chrom, 4)))
}

df_strats <- important_subsets %>%
  map_dfr(~ read_strat(.x[[1]], .x[[2]]))

df_stratified_bed <- df_strats %>%
  left_join(
    df_error_bed,
    join_by(s_chrom == chrom, overlaps(s_start, s_end, start, end))
  )

df_stratified_counts <- df_stratified_bed %>%
  filter(!is.na(error_type)) %>%
  group_by(name, Type) %>%
  summarize(
    nerror = n(),
    .groups = "drop"
  )

df_variants_and_errors <- df_happy %>%
  left_join(df_stratified_counts, by = c(Subset = "name", Type = "Type")) %>%
  mutate(
    FN = phred(TRUTH.FN / Subset.IS_CONF.Size),
    FP = phred(QUERY.FP / Subset.IS_CONF.Size),
    FNFP =phred((TRUTH.FN + QUERY.FP) / Subset.IS_CONF.Size),
    error_pm = phred(nerror / Subset.IS_CONF.Size),
  ) %>%
  filter(Subtype == "*") %>%
  filter(Subset %in% map_chr(important_subsets, ~ .x[[2]])) %>%
  rename(Precision = METRIC.Precision, Recall = METRIC.Recall, F1 = METRIC.F1_Score) %>%
  select(Type, Subset, FP, FN, FNFP, Precision, Recall, F1, error_pm) %>%
  pivot_longer(c(FN, FP, FNFP, Precision, Recall, F1), names_to = "metric", values_to = "metric_value") %>%
  mutate(
    Subset = case_when(
      Subset == "SimpleRepeat_homopolymer_ge12_slop5"     ~ "01 HP (≥12bp)",
      Subset == "SimpleRepeat_homopolymer_ge21_slop5"     ~ "02 HP (≥21bp)",
      Subset == "SimpleRepeat_imperfecthomopolge11_slop5" ~ "03 Imp. HP (≥11bp)",
      Subset == "SimpleRepeat_imperfecthomopolge21_slop5" ~ "04 Imp HP (≥21bp)",
      Subset == "SimpleRepeat_diTR_10to49_slop5"          ~ "05 2-mer (10-49bp)",
      Subset == "SimpleRepeat_triTR_50to149_slop5"        ~ "06 3-mer (50-149bp)",
      Subset == "SimpleRepeat_quadTR_50to149_slop5"       ~ "07 4-mer (50-149bp)",
      Subset == "AllTandemRepeatsandHomopolymers_slop5"   ~ "08 TR/HP",
      Subset == "gclt25orgt65_slop50"                     ~ "09 High/Low GC",
      Subset == "lowmappabilityall"                       ~ "10 Lowmap",
      Subset == "segdups"                                 ~ "11 Segdups",
      Subset == "alldifficultregions"                     ~ "12 All Difficult",
      Subset == "AllAutosomes"                            ~ "13 All Autosomes"
    ) %>%
      factor() %>%
      fct_relabel(~ str_replace(.x, "[0-9]+ ", ""))
  ) %>%
  mutate(Type = if_else(Type == "SNP", "SNP / Subst. Error", "INDEL / INDEL Error"))

# TODO combine INDEL and SNP
df_variants_and_errors %>%
  filter(metric %in% c("FP", "FN", "FNFP")) %>%
  mutate(
    metric = factor(
      if_else(metric == "FNFP", "FN+FP", metric),
      levels = c("FN", "FP", "FN+FP")
    )
  ) %>%
  ggplot(aes(metric_value, error_pm, color = Subset)) +
  geom_abline(slope = 1, intercept = 0, color = "#888888", linetype = "dotted") +
  geom_point() +
  facet_wrap(c("metric", "Type"), nrow = 3) +
  pretty_theme +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 3, byrow = FALSE)) +
  labs(x = "phred(variant label per base)", y = "phred(genome error per base)", color = "Stratification")
ggsave(snakemake@output[["errors_per_base"]], units = "mm", width = 80, height = 140, device = cairo_pdf)

df_variants_and_errors %>%
  filter(metric %in% c("Precision", "Recall", "F1")) %>%
  mutate(metric_value = phred(1 - metric_value)) %>%
  ggplot(aes(metric_value, error_pm, color = Subset)) +
  geom_point() +
  facet_wrap(c("metric", "Type"), nrow = 3) +
  theme(legend.position = "bottom") +
  pretty_theme +
  guides(color = guide_legend(ncol = 3, byrow = FALSE)) +
  labs(x = "phred(1-metric)", y = "phred(genome error per base)", color = "Stratification")
ggsave(snakemake@output[["precision_recall"]], units = "mm", width = 80, height = 140, device = cairo_pdf)

