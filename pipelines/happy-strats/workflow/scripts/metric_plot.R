library(tidyverse)

## TODO wet...
pretty_theme <-
  theme(
    text = element_text(size = 6),
    line = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.25, "mm"),
    axis.line = element_line(linewidth = 0.15),
    legend.box.spacing = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0.75, 0, "mm"),
    plot.margin = margin(1, 1, 1, 1, "mm"),
    legend.key.size = unit(0.55, "lines"),
    legend.spacing = unit(1.5, "mm"),
    strip.text = element_text(size = rel(1), margin = margin(1, 1, 1, 1, "mm")),
    strip.background = element_rect(
      linetype = "blank",
      fill = "gray"
    )
  )

subsets <- c(
  "*" = "All",
  "alldifficultregions" = "Difficult",
  "lowmappabilityall" = "Lowmap",
  "gclt30orgt55_slop50" = "High/Low GC",
  "segdups" = "Segdup",
  "AllTandemRepeats_le50bp_slop5" = "Short TRs",
  "SimpleRepeat_diTR_10to49_slop5" = "2-mer TRs",
  "SimpleRepeat_homopolymer_ge12_GC_slop5" = "HPs (GC)",
  "SimpleRepeat_homopolymer_ge12_AT_slop5" = "HPs (AT)"
  )

read_df <- function(path, which) {
  read_csv(path, col_types = cols(
    Type = "c",
    Subtype = "c",
    Subset = "c",
    Filter = "c",
    METRIC.Recall = "d",
    METRIC.Precision = "d",
    METRIC.F1_Score = "d",
    .default = "-")
    ) %>%
    mutate(which = which) %>%
    rename(Recall = METRIC.Recall,
           Precision = METRIC.Precision,
           F1 = METRIC.F1_Score
           )
}

df <- bind_rows(
  read_df(snakemake@input[["illumina"]], "Illumina"),
  read_df(snakemake@input[["hifi"]], "Hifi")
) %>%
  filter(Filter == "PASS") %>%
  pivot_longer(c(Recall, Precision, F1), names_to = "metric", values_to = "value") %>%
  pivot_wider(id_cols = c(Type, Subtype, Subset, Filter, metric),
              names_from = which,
              values_from = value)

df %>%
  filter(Subtype == "*") %>%
  filter(Subset %in% names(subsets)) %>%
  mutate(Subset = map_chr(Subset, ~ subsets[[.x]])) %>%
  mutate(Subset = factor(Subset, levels = subsets)) %>%
  mutate(Type = if_else(Type == "SNP", "SNV", Type)) %>%
  ggplot(aes(Illumina, Hifi, color = Subset)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point() +
  facet_wrap(c("Type", "metric")) +
  scale_x_continuous(limits = c(0.65, 1)) +
  scale_y_continuous(limits = c(0.65, 1)) +
  pretty_theme
ggsave(snakemake@output[[1]], width = 120, height = 80)
