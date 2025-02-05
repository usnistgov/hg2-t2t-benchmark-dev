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

ensure <- function(test, msg) {
  if (!test) {
    print(test)
    stop(sprintf("ERROR %s", msg))
  }
}

ensure_empty <- function(df, msg) {
  if (nrow(df) > 0) {
    print(df)
    stop(sprintf("ERROR %s", msg))
  }
}

sub_gt <- function(gt, pat) {
  i <- ifelse(pat, 1, 3)
  s <- str_sub(gt, i, i)
  s[s == "."] <- NA
  as.integer(s)
}

# match: perfect match
# collapse: truth het, query hom with one matching allele (1|1 vs 1|2)
# rev_collapse: truth hom, query het with one matching allele (1|2 vs 1|1)
# hom_mismatch: both hom, alleles mismatch (1|1 vs 2|2)
# hom_het_mismatch: truth hom, query het, no matching allele (3|3 vs 1|2)
# het_hom_mismatch: truth het, query hom, no matching allele (2|3 vs 1|1)
# misphase: both het, phase mismatch (1|0 vs 0|1)
# het_mismatch: both het, one allele matches (1|3 vs 1|2)
# misphase_het_mismatch: like het_mismatch but with a flip (1|3 vs 2|1)
# bi_het_mismatch: both het, no matches (3|4 vs 1|2)

gt_to_name <- function(tp, tm, qp, qm) {
  thom <- tp == tm
  qhom <- qp == qm
  exact_matches <- as.integer(tp == qp) + as.integer(tm == qm)
  cross_matches <- as.integer(tp == qm) + as.integer(tm == qp)
  case_when(
    exact_matches == 2 ~ "match",
    thom & qhom & exact_matches == 0 ~ "hom_mismatch",
    thom & !qhom & exact_matches == 1 ~ "rev_collapse",
    thom & !qhom & exact_matches == 0 ~ "hom_het_mismatch",
    !thom & qhom & exact_matches == 1 ~ "collapse",
    !thom & qhom & exact_matches == 0 ~ "het_hom_mismatch",
    !thom & !qhom & exact_matches == 0 & cross_matches == 2 ~ "misphase",
    !thom & !qhom & exact_matches == 1 & cross_matches == 0 ~ "het_mismatch",
    !thom & !qhom & exact_matches == 0 & cross_matches == 1 ~ "misphase_het_mismatch",
    !thom & !qhom & exact_matches == 0 & cross_matches == 0 ~ "bi_het_mismatch",
    TRUE ~ "other"
  )
}

group_indices_ordered <- function(grouped_df) {
  grouped_df %>%
    group_indices() %>%
    match(., unique(.))
}

field_to_matrix <- function(xs, s) {
  ys <- strsplit(xs, s, fixed = TRUE)
  do.call(rbind, lapply(ys, `length<-`, max(lengths(ys))))
}

field_indices <- function(m, f) {
  n <- nrow(m)
  # last column is a dummy to represent NA if f is not found
  i <- apply(cbind(m == f, TRUE), 1, which.max)
  (i - 1) * n + 1:n 
}

select_allele <- function(m, i) {
  n <- nrow(m)
  m[i * n + 1:n]
}

#annotate_vcf <- function(df) {
#  .df <- df %>%
#    filter(!chrom %in% c("chrX", "chrY")) %>%
#    mutate(chrom = as.integer(str_sub(chrom, 4)))
#  .format <- field_to_matrix(.df$format, ":")
#  .truth <- cbind(field_to_matrix(.df$truth, ":"), NA)
#  .query <- cbind(field_to_matrix(.df$query, ":"), NA)
#  .alleles <- cbind(.df$ref, field_to_matrix(.df$alt, ","))
#  gti <- field_indices(.format, "GT")
#  bdi <- field_indices(.format, "BD")
#  bki <- field_indices(.format, "BK")
#  .df %>%
#    mutate(
#      t_gt = .truth[gti],
#      q_gt = .query[gti],
#      t_bd = .truth[bdi],
#      q_bd = .query[bdi],
#      t_bk = .truth[bki],
#      q_bk = .query[bki],
#      tpi = replace_na(sub_gt(t_gt, TRUE), 0),
#      tmi = replace_na(sub_gt(t_gt, FALSE), 0),
#      qpi = replace_na(sub_gt(q_gt, TRUE), 0),
#      qmi = replace_na(sub_gt(q_gt, FALSE), 0),
#      tp = select_allele(.alleles, tpi),
#      tm = select_allele(.alleles, tmi),
#      qp = select_allele(.alleles, qpi),
#      qm = select_allele(.alleles, qmi),
#      variant_class = gt_to_name(
#        replace_na(tpi, 0),
#        replace_na(tmi, 0),
#        replace_na(qpi, 0),
#        replace_na(qmi, 0)
#      )
#    ) %>%
#    select(-format, -truth, -query)
#}

df_hits <- read_tsv(
  snakemake@input[["hit"]],
  na = "NA",
  col_types = cols(
    isec = "i",
    isec_n = "-",
    isec_nv = "i",
    isec_ng = "i",
    chrom = "i",
    vstart = "i",
    vend = "i",
    vid = "i",
    gid = "i",
    adj_realstart = "i",
    adj_realend = "i",
    match_id = "i",
    match_nv = "i",
    match_ng = "i",
    overlap_id = "i",
    genome_expected = "l",
    rowmatch = "l",
    realstart = "i",
    realend = "i",
    .default = "c"
  )
) %>%
  add_column(hits = TRUE) %>%
  # Remove all perfectly matching variants that are part of a match since these
  # might cause some weird stuff later
  select(-vend) %>%
  relocate(hap, chrom, vstart) %>%
  filter(!(
    !is.na(truth_alt) &
      !is.na(query_alt) &
       truth_alt == query_alt &
       !is.na(match_id)
  )) %>%
  group_by(match_id) %>%
  mutate(match_nv1 = length(unique(vid))) %>%
  ungroup()

df_v_nohit <- read_tsv(
  snakemake@input[["v_nohit"]],
  na = "NA",
  col_types = cols(
    chrom = "i",
    vstart = "i",
    vend = "i",
    vid = "i",
    isec = "i",
    isec_n = "-",
    isec_ng = "-",
    isec_nv = "-",
    overlap_id = "-",
    adj_realstart = "-",
    adj_realend = "-",
    genome_expected = "-",
    rowmatch = "l",
    .default = "c"
  )
) %>%
  add_column(nohit = TRUE)

df_g <- read_tsv(
  snakemake@input[["gbench"]],
  na = "NA",
  col_types = cols(
    gchrom = "c",
    gid = "i",
    vid = "i",
    src_hap = "c",
    dst_hap = "c",
    error_type = "c",
    .default = "-"
  )
) %>%
  filter(src_hap == dst_hap) %>%
  rename(hap = src_hap) %>%
  select(-dst_hap)

df_g_nohit <- read_tsv(
  snakemake@input[["g_nohit"]],
  na = "NA",
  col_types = "i--iiic--ci--------"
) %>%
  add_column(hits = FALSE)

proc_df <- function(df) {
  df %>%
    separate_longer_delim(data, "~") %>%
    separate(data, ";", into = c("error_type", NA, NA, NA, NA, NA, "gid", "error_group", "error_group_size", NA)) %>%
    separate(chrom, "_", into = c("chrom", "hap")) %>%
    mutate(gid = as.integer(gid)) %>%
    mutate(error_group_size = as.integer(error_group_size)) %>%
    filter(!chrom %in% c("chrX", "chrY")) %>%
    mutate(chrom = as.integer(str_sub(chrom, 4)))
}

df_errors0 <- read_tsv(
  snakemake@input[["src"]],
  col_names = c("chrom", "start", "end", "data"),
  col_types = "ciic"
) %>%
  proc_df

df_errors1 <- read_tsv(
  snakemake@input[["projectable"]],
  col_names = c("chrom", "start", "end", "overlap", "data"),
  col_types = "ciidc"
) %>%
  proc_df

read_smallvar <- function(which) {
  df_smallvar_pat <- read_tsv(
    snakemake@input[[sprintf("smallvar_%s", which)]],
    col_types = "c--iii----------",
    col_names = c("chrom", "start", "end", "gid")
  ) %>%
  add_column(hap = if_else(which == "pat", "PATERNAL", "MATERNAL"))
}

df_smallvar <- bind_rows(
  read_smallvar("pat"),
  read_smallvar("mat")
)

df_small_ids <- df_smallvar %>%
  select(gid, hap) %>%
  unique()

df_errors1_n <- df_errors1 %>%
  group_by(hap, gid) %>%
  tally()

df_vcf_rid_anno <- read_tsv(
  snakemake@input[["annotated"]],
  col_types = c(
    rid0 = "i",
    rid1 = "i",
    chrom = "i",
    start = "i",
    ref = "c",
    alt = "c",
    tp = "c",
    tm = "c",
    qp = "c",
    qm = "c",
    t_gt = "c",
    q_gt = "c",
    t_bd = "c",
    q_bd = "c",
    t_bk = "c",
    q_bk = "c",
    error_type = "c",
    variant_status = "c",
    variant_class = "c",
    ng_group = "c",
    group_name = "c",
    regions = "c",
    .default = "-"
  )
)

df_errors_grouped <- df_errors0 %>%
  left_join(df_errors1_n, by = c("gid", "hap")) %>%
  mutate(projected = !is.na(n)) %>%
  anti_join(df_small_ids, by = c("gid", "hap")) %>%
  mutate(
    status = case_when(
      projected ~ "Outside Regions",
      !projected ~ "Not Projected",
      TRUE ~ "unknown"
    )
  ) %>%
  replace_na(list(n = 0)) %>%
  select(hap, gid, status, error_type) %>%
  unique()

df_hits %>%
  select(hap, gid, match_id, error_type, hits) %>%
  unique() %>%
  bind_rows(df_g_nohit) %>%
  mutate(
    status = case_when(
      !is.na(match_id) ~ "Exact Match",
      hits ~ "Possible Match",
      TRUE ~ "No Match"
    )
  ) %>%
  replace_na(list(genome_expected = TRUE)) %>%
  bind_rows(df_errors_grouped) %>%
  ggplot(
    aes(
      x = factor(status, levels = c("Exact Match", "Possible Match", "No Match", "Outside Regions", "Not Projected")),
      fill = error_type
    )
  ) +
  geom_bar() +
  theme(legend.position = "bottom") +
  pretty_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = NULL, y = "Number of Genome Errors", fill = "Genome\nError Type")
ggsave(snakemake@output[["g_matches"]], units = "mm", width = 52, height = 60)

df_hit_match_counts <- df_hits %>%
  filter(!is.na(match_id)) %>%
  select(isec, match_id, match_nv, match_ng) %>%
  unique() %>%
  rename(nv = match_nv, ng = match_ng) %>%
  add_column(which = "matched")

df_hit_unmatch_counts <- df_hits %>%
  filter(is.na(match_id)) %>%
  select(isec, vid, gid) %>%
  group_by(isec) %>%
  summarize(
    nv = length(unique(vid)),
    ng = length(unique(gid))
  ) %>%
  add_column(which = "unmatched")

df_hit_counts <- bind_rows(
  df_hit_match_counts,
  df_hit_unmatch_counts
) %>%
  arrange(isec, match_id) %>%
  # NOTE: id's are not necessarily unique globally but will be unique within
  # the match/unmatch groups
  mutate(id = coalesce(match_id, isec))

# lots of variants (that actually match) only have one genome error
df_hit_counts %>%
  mutate(which = if_else(which == "matched", "Exact Match", "Possible Match")) %>%
  ggplot(aes(ng, fill = factor(nv))) +
  geom_bar() +
  facet_wrap(c("which")) +
  pretty_theme +
  labs(
    # TODO this is actually a bit misleading since the matches aren't per intersection
    x = "genome errors\nper intersection",
    y = "intersection count",
    fill = "variants per\nintersection"
  )
ggsave(snakemake@output[["vg_matches"]], units = "mm", width = 80, height = 80)

df_vcf_rid_anno %>%
  filter(!is.na(variant_status)) %>%
  # missing_hit is so small that they might as well be "no match"
  mutate(variant_status = case_when(
    variant_status == "no_hit" ~ "No Match",
    variant_status == "missing_hit" ~ "No Match",
    variant_status == "matched" ~ "Exact Match",
    variant_status == "unmatched" ~ "Possible Match",
    TRUE ~ "fixmerightmeow"
  )) %>%
  mutate(variant_status = fct_relevel(factor(variant_status), "No Match", after = 2)) %>%
  ggplot(aes(variant_status, fill = ng_group)) +
  geom_bar() +
  labs(x = NULL, y = "Number of Variants", fill = "Genome Error\nCount (pat, mat)") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  pretty_theme
ggsave(snakemake@output[["v_matches"]], units = "mm", width = 55, height = 60)

df_vcf_rid_anno %>%
  filter(!is.na(rid0)) %>%
  ggplot(aes(y = variant_class, fill = error_type)) +
  geom_bar() +
  facet_wrap(c("group_name")) +
  pretty_theme +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = NULL, fill = "Genome\nError Type")
ggsave(snakemake@output[["class_by_type"]], units = "mm", width = 90, height = 130)

df_vcf_rid_anno %>%
  filter(!is.na(rid0)) %>%
  ggplot(aes(y = variant_class, fill = ng_group)) +
  geom_bar() +
  facet_wrap(c("group_name")) +
  pretty_theme +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = NULL, fill = "Genome Error\nCount (pat, mat)")
ggsave(snakemake@output[["class_by_count"]], units = "mm", width = 90, height = 130)

# collapse variant classes into "umbrella classes" which better describe the
# "mechanism" for how the variant arose
df_vcf_rid_anno_short <- df_vcf_rid_anno %>%
  filter(!(variant_class == "match" & group_name == "unmatched")) %>%
  filter(!is.na(rid0)) %>%
  mutate(
    variant_class = case_when(
      variant_class == "misphase" ~ "Misphases",
      variant_class == "other" ~ "Other Errors",
      variant_class == "misphase_het_mismatch" ~ "Misphases/Seq Errors",
      variant_class %in% c("collapse", "het_hom_mismatch",
                           "collapse/het_hom_mismatch") ~ "Collapses",
      variant_class %in% c("rev_collapse", "hom_het_mismatch", "het_mismatch",
                           "bi_het_mismatch", "hom_mismatch",
                           "bi_het_mismatch/het_mismatch") ~ "Seq Errors",
    )
  ) %>%
  mutate(variant_class = factor(variant_class, levels = c("Collapses", "Seq Errors", "Misphases", "Misphases/Seq Errors", "Other Errors")))

df_vcf_rid_anno_short %>%
  filter(group_name != "unmatched") %>%
  mutate(error_type = fct_relevel(factor(error_type), "MIXED", after = 2)) %>%
  ggplot(aes(y = fct_rev(variant_class), fill = error_type)) +
  geom_bar() +
  facet_wrap(c("group_name")) +
  labs(x = NULL, y = NULL, fill = "Genome Error\nType") +
  theme(legend.position = "bottom") +
  pretty_theme
ggsave(snakemake@output[["class_by_type_short"]], width = 80, height = 57, units = "mm")

df_vcf_rid_anno_short %>%
  filter(group_name != "unmatched") %>%
  ggplot(aes(y = fct_rev(variant_class), fill = ng_group)) +
  geom_bar() +
  facet_wrap(c("group_name")) +
  labs(x = NULL, y = NULL, fill = "Genome Error\nCount (pat, mat)") +
  theme(legend.position = "bottom") +
  pretty_theme
ggsave(snakemake@output[["class_by_count_short"]], width = 80, height = 57, units = "mm")

#
# how many het_hom_mismatches are "average collapses"? (needed to test the
# assumption we use above that allows us to group het_hom_mismatches with
# collapses)
#

df_ave_collapse1 <- df_vcf_rid_anno %>%
  filter(variant_class == "het_hom_mismatch") %>%
  filter(!group_name %in% c("unmatched", "complex")) %>%
  # just look at 1 line for now
  filter(str_detect(group_name, "1 Line")) %>%
  mutate(
    tplen = str_length(tp),
    tmlen = str_length(tm),
    tshorter = if_else(tplen < tmlen, tplen, tmlen),
    tlonger = if_else(tplen < tmlen, tmlen, tplen),
    # ASSUME both are the same
    qlen = str_length(qp),
    is_ave_collapse = tshorter < qlen & qlen < tlonger
  )

df_ave_collapse2 <- df_vcf_rid_anno %>%
  filter(variant_class == "het_hom_mismatch") %>%
  filter(!group_name %in% c("unmatched", "complex")) %>%
  filter(str_detect(group_name, "2 Line")) %>%
  arrange(rid1, t_gt) %>%
  ## select(-matches("_(bd|bk)"), -ends_with("i"))
  group_by(rid1, group_name) %>%
  summarize(
    chrom = first(chrom),
    start = first(start),
    ref0 = first(ref),
    # ASSUME both are the same
    q0 = first(qp),
    ref1 = last(ref),
    tp1 = last(tp),
    tm1 = last(tm),
    .groups = "drop"
  ) %>%
  select(-rid1) %>%
  mutate(
    qindel = str_length(q0) - str_length(ref0),
    tp_indel = str_length(tp1) - str_length(ref1),
    tm_indel = str_length(tm1) - str_length(ref1),
    tshorter = if_else(tp_indel < tm_indel, tp_indel, tm_indel),
    tlonger = if_else(tp_indel < tm_indel, tm_indel, tp_indel),
    is_ave_collapse = tshorter < qindel & qindel < tlonger
  )

df_ave_collapse1 %>%
  select(group_name, is_ave_collapse) %>%
  bind_rows(df_ave_collapse2 %>% select(group_name, is_ave_collapse)) %>%
  mutate(y = if_else(is_ave_collapse, "Ave. Collapse", "Not Ave. Collapse")) %>%
  ggplot(aes(y = y, fill = group_name)) +
  geom_bar() +
  pretty_theme +
  labs(x = "Number of Variants", y = NULL, fill = "category")
ggsave(snakemake@output[["ave_collapse"]], units = "mm", width = 80, height = 80)

filter_regions <- function(xs) {
  redundant <- xs %in% c("AllAutosomes", "CONF", "TS_contained")
  genome <- str_detect(xs, "^HG00")
  ancestry <- str_detect(xs, "^ancestry_")
  xs[!(redundant | genome | ancestry)]
}

df_seq <- df_vcf_rid_anno %>%
  mutate(
    big_variant_class = case_when(
      variant_class %in% c("hom_mismatch", "het_hom_mismatch",
                           "bi_het_mismatch/het_mismatch", "bi_het_mismatch",
                           "het_mismatch", "rev_collapse", "hom_het_mismatch") ~
        "Seq Errors",
      variant_class == "misphase" ~ "Misphases",
      variant_class == "other" ~ "Other Errors",
      variant_class == "misphase_het_mismatch" ~ "Seq Errors/Misphases",
      variant_class %in% c("collapse", "het_hom_mismatch",
                           "collapse/het_hom_mismatch") ~
        "Collapses",
      TRUE ~ NA
    )
  ) %>%
  filter(!is.na(big_variant_class)) %>%
  select(chrom, start, ref, alt, tp, tm, qp, qm, variant_class,
         big_variant_class, variant_status, ng_group, regions) %>%
  mutate(
    tplen = str_length(tp),
    tmlen = str_length(tm),
    qplen = str_length(qp),
    qmlen = str_length(qm),
    pdiff = qplen - tplen,
    mdiff = qmlen - tmlen,
    # ASSUME no SVs
    ptype = case_when(
      tp == qp ~ NA,
      tplen == 1 & qplen == 1 ~ "Subst. Error",
      TRUE ~ "INDEL Error"
    ),
    mtype = case_when(
      tm == qm ~ NA,
      tmlen == 1 & qmlen == 1 ~ "Subst. Error",
      TRUE ~ "INDEL Error"
    ),
    type = case_when(
      is.na(ptype) & is.na(mtype) ~ NA,
      is.na(ptype) ~ mtype,
      is.na(mtype) ~ ptype,
      ptype == mtype ~ ptype,
      TRUE ~ "MIXED"
    )
  ) %>%
  mutate(regions = str_replace(regions, "BS=[0-9]+;Regions=", "")) %>%
  mutate(regions = map(str_split(regions, ","), filter_regions)) %>%
  mutate(
    inHP = map_lgl(regions, ~ any(str_detect(.x, "^SimpleRepeat_homopolymer_"))),
    #inHPshort = map_lgl(regions, ~ any(str_detect(.x, "SimpleRepeat_homopolymer_4to6"))),
    #inHPmid = map_lgl(regions, ~ any(str_detect(.x, "SimpleRepeat_homopolymer_7to11"))),
    #inHPlong = map_lgl(regions, ~ any(str_detect(.x, "SimpleRepeat_homopolymer_ge12"))),
    #inHPxlong = map_lgl(regions, ~ any(str_detect(.x, "SimpleRepeat_homopolymer_ge21"))),
    inHPTR = map_lgl(regions, ~ any(str_detect(.x, "^AllTandemRepeatsandHomopolymers"))),
    #inGC = map_lgl(regions, ~ any(str_detect(.x, "^gc")))
  ) %>%
  mutate(
    region_class = case_when(
      inHPTR & !inHP ~ "TRs",
      inHP ~ "HPs",
      TRUE ~ "Other"
    )
  )

df_seq %>%
  mutate(region_class = fct_relevel(region_class, "Other", after = 2) %>% fct_rev()) %>%
  mutate(big_variant_class = fct_relevel(big_variant_class, c("Collapses", "Seq Errors", "Misphases", "Seq Errors/Misphases"))) %>%
  mutate(type = fct_relevel(type, "MIXED", after = 2)) %>%
  ggplot(aes(y = region_class, fill = type)) +
  geom_bar() +
  facet_wrap(c("big_variant_class"), ncol = 2) +
  pretty_theme +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Number of Variants", y = NULL, fill = "Genome Error\nModality")
ggsave(snakemake@output[["errors_by_region"]], units = "mm", width = 60, height = 60)

