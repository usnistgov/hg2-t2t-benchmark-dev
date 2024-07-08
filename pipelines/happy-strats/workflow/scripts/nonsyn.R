library(tidyverse)
library(ggpubr)

read_vcf <- function(path) {
  vcf_cols <- c(
    "chrom",
    "pos",
    "id",
    "ref",
    "alt",
    "qual",
    "filter",
    "info",
    "keys",
    "values"
  )
  vcf_types <- "cicccicccc"
  read_tsv(path, col_names = vcf_cols, col_types = vcf_types) %>%
    mutate(
      ref_length = str_length(ref),
      alt_length = str_length(alt),
      indel_length = alt_length - ref_length,
      vartype = case_when(
        ref_length == 1 & alt_length == 1 ~ "SNV",
        abs(indel_length) < 50 ~ "INDEL",
        abs(indel_length) >= 50 ~ "SV",
        TRUE ~ NA,
        ),
      chridx = case_when(
        chrom == "chrX" ~ 23,
        chrom == "chrY" ~ 24,
        TRUE ~ as.integer(str_extract(chrom, "chr([0-9]+)", 1)),
        ),
      ) %>%
    separate_wider_delim(values, delim = ":", names = c("GT", "AD")) %>%
    select(-keys, -id, -filter, -qual)
}

vcf_cols <- c(
  "chrom",
  "pos",
  "id",
  "ref",
  "alt",
  "qual",
  "filter",
  "info",
  "keys",
  "values"
)

vcf_types <- "cicccicccc"

syn <- read_vcf(snakemake@input[["syn"]])

nonsyn <- read_vcf(snakemake@input[["nonsyn"]])

df <- bind_rows(
  mutate(syn, origin = "syn"),
  mutate(nonsyn, origin = "nonsyn")
) %>%
  arrange(chridx, pos, ref) %>%
  filter(chridx != 24)

nvars <- df %>%
  filter(origin == "nonsyn") %>%
  group_by(chrom, chridx, vartype) %>%
  tally() %>%
  ggplot(
    aes(
      n,
      fct_rev(fct_reorder(chrom, chridx)),
      fill = factor(vartype, levels = c("SV", "INDEL", "SNV"))
    )
  ) +
  geom_col() +
  labs(
    x = "Number of NonSyn Variants",
    y = NULL,
    fill = "Variant Type"
  )

pvars <- df %>%
  group_by(chrom, chridx, vartype, origin) %>%
  tally() %>%
  pivot_wider(c(chrom, chridx, vartype), names_from = origin, values_from = n) %>%
  mutate(frac_nonsyn = nonsyn / (syn + nonsyn)) %>%
  ggplot(
    aes(
      frac_nonsyn,
      fct_rev(fct_reorder(chrom, chridx)),
      fill = factor(vartype, levels = c("SV", "INDEL", "SNV")),
    )
  ) +
  geom_col(position = "dodge") +
  scale_x_continuous(labels = scales::percent) +
  ## facet_wrap(c("vartype"), scales = "free") +
  labs(
    x = "Perc. Increase N Variants\n(NonSyn / All)",
    y = NULL,
    fill = "Variant Type"
  ) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

nbases <- df %>%
  filter(origin == "nonsyn") %>%
  ## mutate(vartype = if_else(vartype == "SV", "SV", "SNV+INDEL")) %>%
  group_by(chrom, chridx, vartype) %>%
  summarize(
    total_bases = sum(ref_length)
  ) %>%
  ggplot(
    aes(
      total_bases,
      fct_rev(fct_reorder(chrom, chridx)),
      fill = factor(vartype, levels = c("SV", "INDEL", "SNV")),
    )
  ) +
  geom_col(position = "dodge") +
  scale_x_log10() +
  ## facet_wrap(c("vartype"), scales = "free") +
  labs(
    x = "Number of NonSyn\nVariant Bases",
    y = NULL,
    fill = "Variant Type",
    ) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggarrange(nvars, pvars, nbases, common.legend = TRUE, nrow = 1, align = "h",
          widths = c(1, 0.9, 0.9))
ggsave(snakemake@output[[1]])
