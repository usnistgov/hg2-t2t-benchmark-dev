library(tidyverse)

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

annotate_vcf <- function(df) {
  .df <- df %>%
    filter(!chrom %in% c("chrX", "chrY")) %>%
    mutate(chrom = as.integer(str_sub(chrom, 4)))
  .format <- field_to_matrix(.df$format, ":")
  .truth <- cbind(field_to_matrix(.df$truth, ":"), NA)
  .query <- cbind(field_to_matrix(.df$query, ":"), NA)
  .alleles <- cbind(.df$ref, field_to_matrix(.df$alt, ","))
  gti <- field_indices(.format, "GT")
  bdi <- field_indices(.format, "BD")
  bki <- field_indices(.format, "BK")
  .df %>%
    mutate(
      t_gt = .truth[gti],
      q_gt = .query[gti],
      t_bd = .truth[bdi],
      q_bd = .query[bdi],
      t_bk = .truth[bki],
      q_bk = .query[bki],
      tpi = replace_na(sub_gt(t_gt, TRUE), 0),
      tmi = replace_na(sub_gt(t_gt, FALSE), 0),
      qpi = replace_na(sub_gt(q_gt, TRUE), 0),
      qmi = replace_na(sub_gt(q_gt, FALSE), 0),
      tp = select_allele(.alleles, tpi),
      tm = select_allele(.alleles, tmi),
      qp = select_allele(.alleles, qpi),
      qm = select_allele(.alleles, qmi),
      variant_class = gt_to_name(
        replace_na(tpi, 0),
        replace_na(tmi, 0),
        replace_na(qpi, 0),
        replace_na(qmi, 0)
      )
    ) %>%
    select(-format, -truth, -query)
}

df_vcf0 <- read_tsv(
  snakemake@input[["vcf"]],
  col_types = "ciicc---ccc",
  col_names = c("chrom", "start", "vid", "ref", "alt", "format", "truth", "query")
)

df_vcf <- df_vcf0 %>%
  annotate_vcf()

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

#
# some overall plots
#

df_hits %>%
  select(hap, gid, match_id, error_type, hits) %>%
  unique() %>%
  bind_rows(df_g_nohit) %>%
  mutate(
    status = case_when(
      !is.na(match_id) ~ "matched",
      hits ~ "unmatched",
      TRUE ~ "no_hit"
    )
  ) %>%
  replace_na(list(genome_expected = TRUE)) %>%
  ggplot(
    aes(
      x = factor(status, levels = c("matched", "unmatched", "missing_hit", "no_hit")),
      fill = error_type
    )
  ) +
  geom_bar() +
  labs(x = NULL, y = "number of genome errors", fill = "expected variant hit")
ggsave(snakemake@output[["g_matches"]])

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
  ggplot(aes(ng, fill = factor(nv))) +
  geom_bar() +
  facet_wrap(c("which")) +
  labs(
    # TODO this is actually a bit misleading since the matches aren't per intersection
    fill = "variants per\nintersection",
    y = "intersection count",
    x = "genome errors\nper intersection"
  )
ggsave(snakemake@output[["vg_matches"]])

#
# (attempt to) classify stuff
#

df_error_types <- df_hits %>%
  mutate(.is_consensus = error_type == "CONSENSUS") %>%
  group_by(match_id, isec) %>%
  summarize(
    error_type = case_when(
      all(.is_consensus) ~ "CONSENSUS",
      all(!.is_consensus) ~ "PHASING",
      TRUE ~ "MIXED"
    ),
    .groups = "drop"
  )

# for each vid, pivot all genome errors onto pat/mat haplotypes; each row is
# one vid
df_isec_match <- df_hits %>%
  select(hap, vid, match_id, match_nv, match_nv1, match_ng, isec, isec_ng) %>%
  unique() %>%
  pivot_wider(
    id_cols = vid,
    values_from = c(match_id, match_nv, match_nv1, match_ng, isec, isec_ng),
    names_from = hap,
    names_glue = "{hap}_{.value}"
  ) %>%
  left_join(
    rename_with(df_error_types, ~ sprintf("pat_%s", .x)),
    by = c("pat_match_id", "pat_isec")
  ) %>%
  left_join(
    rename_with(df_error_types, ~ sprintf("mat_%s", .x)),
    by = c("mat_match_id", "mat_isec")
  ) %>%
  rename_with(~ str_replace(.x, "_match", ""), matches("(nv|ng)$")) %>%
  rename_with(~ str_replace(.x, "match_", "m"), ends_with("match_id")) %>%
  mutate(
    in_match = !is.na(mat_mid) | !is.na(pat_mid),
    error_type = case_when(
      is.na(mat_error_type) ~ pat_error_type,
      is.na(pat_error_type) ~ mat_error_type,
      pat_error_type != mat_error_type ~ "MIXED",
      pat_error_type == mat_error_type ~ pat_error_type,
      TRUE ~ "SOMETHING"
    )
  ) %>%
  replace_na(
    list(
      pat_nv = 0, mat_nv = 0,
      pat_ng = 0, mat_ng = 0,
      pat_isec_ng = 0, mat_isec_ng = 0
    )
  )

# variants without a hit, one row = one vid
df_v_nohit_unique <- df_v_nohit %>%
  mutate(missing_hit = !is.na(isec)) %>%
  select(vid, missing_hit) %>%
  unique() %>%
  add_column(nohit = TRUE)

error_num <- function(p, m) {
  case_when(
    is.na(p) & is.na(m) ~ "none",
    p == 0 & m == 0 ~ "none",
    p > 1 | m > 1 ~ "many",
    TRUE ~ sprintf("%s,%s", as.character(p), as.character(m))
  )
}

# vcf file annotated with genome error match data (or lack thereof); one
# line = one variant
df_vcf_isec_match <- df_vcf %>%
  left_join(df_isec_match, by = "vid") %>%
  left_join(df_v_nohit_unique, by = "vid") %>%
  replace_na(list(nohit = FALSE)) %>%
  mutate(
    ng_group = case_when(
      in_match ~ error_num(pat_ng, mat_ng),
      !in_match ~ error_num(pat_isec_ng, mat_isec_ng),
      TRUE ~ "none"
    ),
    group_name = sprintf(
      "1 Line - %s | %s",
      if_else(t_bd == ".", "none", t_bd),
      if_else(q_bd == ".", "none", q_bd)
    ),
    has_error = !(t_bd == "TP" | q_bd == "TP")
  )

select_pat_or_mat <- function(df, hap0, hap1, i) {
  df_i <- df %>%
    select({{ hap0 }}, {{ i }}) %>%
    filter(!is.na({{ hap0 }})) %>%
    unique() %>%
    rename({{ hap1 }} := {{ i }})
  # each pat/mat id should only appear once
  df_i %>%
    group_by({{ hap0 }}) %>%
    tally() %>%
    filter(n > 1) %>%
    ensure_empty("hap must appear only once")
  df_i
}

# Figure out which indices across pat/mat columns "belong together" (ie have
# shared indices) and index each of the resulting groups. For example, [NA, 1,
# 2] and [1, 2, NA] in three rows across pat/mat columns respectively would be
# grouped together.
group_indices <- function(df, pat, mat, i, c) {
  .pat <- rlang::as_name(enquo(pat))
  .mat <- rlang::as_name(enquo(mat))
  df_mp <- df %>%
    select({{ pat }}, {{ mat }})
  df_pairs <- df_mp %>%
    filter(!is.na({{ pat }}) & !is.na({{ mat }})) %>%
    unique()
  # ensure each of these columns are monotonically increasing, otherwise the
  # indexing thing we do won't work
  df_pairs %>%
    filter(lag({{ pat }}) > {{ pat }}) %>%
    ensure_empty("pat must be in order")
  df_pairs %>%
    filter(lag({{ mat }}) > {{ mat }}) %>%
    ensure_empty("mat must be in order")
  df_i <- df_pairs %>%
    mutate(
      {{ i }} := cumsum(
        lag({{ pat }}, default = 0) < {{ pat }} &
          lag({{ mat }}, default = 0) < {{ mat }}
      )
    )
  df_pat_only <- df_mp %>%
    anti_join(df_pairs, by = .pat) %>%
    filter(!is.na({{ pat }})) %>%
    unique()
  df_mat_only <- df_mp %>%
    anti_join(df_pairs, by = .mat) %>%
    filter(!is.na({{ mat }})) %>%
    unique()
  max_i <- df_i %>% pull({{ i }}) %>% max()
  df_i_all <- bind_rows(df_pat_only, df_mat_only) %>%
    mutate({{ i }} := row_number() + max_i) %>%
    bind_rows(df_i)
  df_i_p <- select_pat_or_mat(df_i_all, {{ pat }}, .pat_i, {{ i }})
  df_i_m <- select_pat_or_mat(df_i_all, {{ mat }}, .mat_i, {{ i }})
  df %>%
    left_join(df_i_p, by = .pat) %>%
    left_join(df_i_m, by = .mat) %>%
    # ASSUME that these will always be the same if both are present, the reason
    # why this is necessary is to fill the rows where only pat or mat was
    # present
    mutate({{ i }} := coalesce(.pat_i, .mat_i)) %>%
    select(-.mat_i, -.pat_i) %>%
    group_by({{ i }}) %>%
    mutate({{ c }} := if_else(is.na({{ i }}), NA, n())) %>%
    ungroup()
}

# variants that are grouped based on the errors they intersect. One line = one
# variant. rid0 is the group based on the intersection index (ie a "hit") and
# rid1 is the same for shared match_ids
df_vcf_rid <- df_vcf_isec_match %>%
  group_indices(pat_isec, mat_isec, rid0, nr0) %>%
  group_indices(pat_mid, mat_mid, rid1, nr1) %>%
  group_by(rid1) %>%
  ungroup() %>%
  mutate(
    nr1 = if_else(in_match, nr1, 0),
  ) %>%
  relocate(vid, rid0, rid1, nr0, nr1)

# Groups of variants with two variants per group, with updated variant classes
# based on the genotypes (and other stuff) of each variant in the pair. The
# general strategy is to imagine what the genotype would look like if the
# variant were written on one line, which in many cases is easy to figure out or
# assume. Note, we don't care about anything beyond paired variants (or single
# which don't require special treatment) since these are insanely complex and
# not that common.
df_ref_matches2 <- df_vcf_rid %>%
  filter(nr1 == 2) %>%
  filter(pat_ng < 2 & mat_ng < 2) %>%
  group_by(rid1) %>%
  # All of these are "split" across 2 variants, so arrange the variants using
  # their BD field so that the first and second are predictable, and therefore
  # we can classify based on their patterns
  arrange(t_bd, q_bd) %>%
  summarize(
    chrom = first(chrom),
    start0 = first(start),
    start1 = last(start),
    ref0 = first(ref),
    ref1 = last(ref),
    tgt0 = first(t_gt),
    tgt1 = last(t_gt),
    qgt0 = first(q_gt),
    qgt1 = last(q_gt),
    tbk0 = first(t_bk),
    tbk1 = last(t_bk),
    qbk0 = first(q_bk),
    qbk1 = last(q_bk),
    tp0 = first(tp),
    tp1 = last(tp),
    tm0 = first(tm),
    tm1 = last(tm),
    qp0 = first(qp),
    qp1 = last(qp),
    qm0 = first(qm),
    qm1 = last(qm),
    nFN = sum(t_bd == "FN"),
    nFP = sum(q_bd == "FP"),
    ntTP = sum(t_bd == "TP"),
    nqTP = sum(t_bd == "TP"),
    pat_ng0 = first(pat_ng),
    mat_ng0 = first(mat_ng),
    pat_ng1 = last(pat_ng),
    mat_ng1 = last(mat_ng)
  ) %>%
  mutate(
    .variant_class = case_when(
      # ASSUME these are collapses where both lines are far away from each other
      # and if written on one line would be 1|2+1|1
      tgt0 == "0|1" & tgt1 == "1|0" & qgt0 == "." & qgt1 == "1|1" ~ "collapse",
      tgt0 == "1|0" & tgt1 == "0|1" & qgt0 == "." & qgt1 == "1|1" ~ "collapse",
      # ASSUME this is basically the above case but written in a strange way
      tgt0 == ".|1" & tgt1 == "1|0" & qgt0 == "." & qgt1 == "1|1" ~ "collapse",
      tgt0 == "1|." & tgt1 == "0|1" & qgt0 == "." & qgt1 == "1|1" ~ "collapse",
      # ASSUME these would be written as 2|1+3|0 on one line (I suppose there
      # is a small chance that it may be valid as 2|1+1|0); if these are bi_het
      # then they should also have two errors on either hap. The alternatives
      # seem to be rare and complicated, either involving other variant lines
      # or genome errors that were missed
      tgt0 == "." & tgt1 == "2|1" & qgt0 == "1|0" & qgt1 == "." ~
        if_else(pat_ng1 == 1 & mat_ng1 == 1, "bi_het_mismatch", "TODO1"),
      tgt0 == "." & tgt1 == "1|2" & qgt0 == "0|1" & qgt1 == "." ~
        if_else(pat_ng1 == 1 & mat_ng1 == 1, "bi_het_mismatch", "TODO1"),
      # ASSUME if these groups intersect with 2 variants then they are bi-het
      # mismatches and the GT fields are like 2|1+0|3 (or reverse). If they
      # only have one error, then either they are also bi_het_mismatches and
      # the error didn't match or liftover or they are just het_mismatches
      # and are actually on 3 lines with one line missing, so either way they
      # don't belong to the "simple" groups
      tgt0 == "." & tgt1 == "2|1" & qgt0 == "0|1" & qgt1 == "." ~
        if_else(pat_ng1 == 1 & mat_ng1 == 1, "bi_het_mismatch", "TODO2"),
      tgt0 == "." & tgt1 == "1|2" & qgt0 == "1|0" & qgt1 == "." ~
        if_else(pat_ng1 == 1 & mat_ng1 == 1, "bi_het_mismatch", "TODO2"),
      # ASSUME these would be written as 1|1+1|2 on single line
      tgt0 == "." & tgt1 == "1|1" & qgt0 == "0|1" & qgt1 == "1|0" ~ "rev_collapse",
      tgt0 == "." & tgt1 == "1|1" & qgt0 == "1|0" & qgt1 == "0|1" ~ "rev_collapse",
      # These are either het_mismatches or bi_het_mismatches, not clear what
      # distinguishes them, they are quite rare regardless
      tgt0 == "." & tgt1 == ".|1" & qgt0 == "0|1" & qgt1 == "." ~ "TODO3",
      tgt0 == "." & tgt1 == "1|." & qgt0 == "1|0" & qgt1 == "." ~ "TODO3",
      # ASSUME these would be written like 1|2+1|1; most of these seem to be
      # split because the half call variant represents a het->hom SNP and the
      # het->null variant is a single base deletion (TODO are these special?)
      tgt0 == "1|0" & tgt1 == ".|1" & qgt0 == "." & qgt1 == "1|1" ~ "collapse",
      tgt0 == "0|1" & tgt1 == "1|." & qgt0 == "." & qgt1 == "1|1" ~ "collapse",
      # obvious case, should be 1|1+2|2 if on the same line (and since they
      # aren't they clearly don't match, hence the different allele numbers)
      tgt0 == "." & tgt1 == "1|1" & qgt0 == "1|1" & qgt1 == "." ~ "hom_mismatch",
      # ASSUME if these were written on the same line that they would look like
      # 1|2+1|3 (ie one het mismatches). This is justified because if the 2 and
      # 3 alleles matched, then the lines wouldn't be split apart like this
      tgt0 == "." & tgt1 == "1|2" & qgt0 == "0|1" & qgt1 == "1|0" ~ "het_mismatch",
      tgt0 == "." & tgt1 == "2|1" & qgt0 == "1|0" & qgt1 == "0|1" ~ "het_mismatch",
      # ASSUME same as above, these would be written on one line like 0|1+0|2
      tgt0 == "." & tgt1 == "0|1" & qgt0 == "0|1" & qgt1 == "." ~ "het_mismatch",
      tgt0 == "." & tgt1 == "1|0" & qgt0 == "1|0" & qgt1 == "." ~ "het_mismatch",
      # Simply comparing the GT field is not enough to classify these, since if
      # compressed to one line they could be written as 1|2+1|1, 1|2+2|2, or
      # 1|2+3|3 (the first of which are collapses). However, many (but not all
      # of these) of these variant pairs start/end on the same position, so we
      # can directly compare the alleles in the vcf to figure this out.
      tgt0 == "." & tgt1 %in% c("1|2", "2|1") & qgt0 == "1|1" & qgt1 == "." ~
        "collapse/het_hom_mismatch",
      # ASSUME that the 1|1 would be written like 2|2 if these were on the same line
      tgt0 == "." & tgt1 %in% c("1|0", "0|1") & qgt0 == "1|1" & qgt1 == "." ~ "het_hom_mismatch",
      TRUE ~ "other"
    )
  ) %>%
  mutate(
    .truth_name = case_when(
      ntTP == 0 & nFN == 0 ~ "none",
      ntTP == 0 & nFN == 1 ~ "FN",
      ntTP == 1 & nFN == 0 ~ "TP",
      ntTP == 1 & nFN == 1 ~ "TP+FN",
      ntTP == 2 & nFN == 0 ~ "2xTP",
      ntTP == 0 & nFN == 2 ~ "2xFN",
      TRUE ~ "other"
    ),
    .query_name = case_when(
      nqTP == 0 & nFP == 0 ~ "none",
      nqTP == 0 & nFP == 1 ~ "FP",
      nqTP == 1 & nFP == 0 ~ "TP",
      nqTP == 1 & nFP == 1 ~ "TP+FP",
      nqTP == 2 & nFP == 0 ~ "2xTP",
      nqTP == 0 & nFP == 2 ~ "2xFP",
      TRUE ~ "other"
    ),
    .group_name = sprintf("2 Line - %s | %s", .truth_name, .query_name)
  )

# update vcf annotations with new variant classes from above using paired
# variants
df_vcf_rid_anno <- df_ref_matches2 %>%
  select(rid1, .variant_class, .group_name) %>%
  right_join(df_vcf_rid, by = "rid1") %>%
  mutate(
    variant_class = coalesce(.variant_class, variant_class),
    group_name = coalesce(.group_name, group_name),
  ) %>%
  select(-starts_with(".")) %>%
  # these are actually het_mismatches but are misclassified as het_hom_mismatches
  mutate(
    .correct_het_hom = group_name == "1 Line - FN | none" &
      ng_group %in% c("1,0", "0,1") &
      variant_class == "het_hom_mismatch",
    variant_class = if_else(.correct_het_hom, "het_mismatch", variant_class),
    group_name = if_else(.correct_het_hom, "2 Line - FN | FP", group_name)
  ) %>%
  mutate(
    variant_class = if_else(nr1 > 2, "other", variant_class),
    group_name = case_when(
      nr1 == 0 ~ "unmatched",
      nr1 > 2 ~ "complex",
      TRUE ~ group_name
    ),
    group_name = case_when(
      variant_class %in% c("TODO1", "TODO2", "TODO3") ~ "complex",
      group_name %in% c(
        "1 Line - none | N",
        "1 Line - none | TP",
        "1 Line - TP | FP",
        "1 Line - TP | none",
        "2 Line - 2xTP | 2xTP",
        "2 Line - FN | none",
        "2 Line - none | none",
        "2 Line - TP | TP",
        "2 Line - TP | TP+FP",
        "2 Line - TP+FN | TP",
        "2 Line - TP+FN | TP+FP"
      ) ~ "complex",
      TRUE ~ group_name
    ),
    variant_class = case_when(
      str_detect(variant_class, "TODO") ~ "bi_het_mismatch/het_mismatch",
      TRUE ~ variant_class
    ),
    variant_status = factor(
      case_when(
        nr1 > 0 ~ "matched",
        nr1 == 0 ~ "unmatched",
        is.na(nr1) & nohit & missing_hit ~ "missing_hit",
        is.na(nr1) & nohit & !missing_hit ~ "no_hit",
        TRUE ~ NA
      ),
      levels = c("matched", "unmatched", "missing_hit", "no_hit")
    )
  )

# print annotated variants for subsequent viewing pleasure
df_vcf_rid_anno %>%
  filter(!is.na(variant_status)) %>%
  select(chrom, start, ref, alt, t_gt, q_gt, t_bd, q_bd, pat_error_type, mat_error_type, error_type, ng_group, variant_class, group_name) %>%
  write_tsv(snakemake@output[["annotated"]])


#
# plot the universe
#

df_vcf_rid_anno %>%
  filter(!is.na(variant_status)) %>%
  ggplot(aes(variant_status, fill = ng_group)) +
  geom_bar() +
  labs(x = NULL, y = NULL, fill = "Error Count\n(pat, mat)")
ggsave(snakemake@output[["v_matches"]])

df_vcf_rid_anno %>%
  filter(!is.na(rid0)) %>%
  ggplot(aes(y = variant_class, fill = error_type)) +
  geom_bar() +
  facet_wrap(c("group_name")) +
  labs(x = NULL, y = NULL, fill = "Error Type")
ggsave(snakemake@output[["class_by_type"]], width = 11, height = 10)

df_vcf_rid_anno %>%
  filter(!is.na(rid0)) %>%
  ggplot(aes(y = variant_class, fill = ng_group)) +
  geom_bar() +
  facet_wrap(c("group_name")) +
  labs(x = NULL, y = NULL, fill = "Error count\n(pat, mat)")
ggsave(snakemake@output[["class_by_count"]], width = 11, height = 10)

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
  ggplot(aes(y = fct_rev(variant_class), fill = error_type)) +
  geom_bar() +
  facet_wrap(c("group_name")) +
  labs(x = NULL, y = NULL, fill = "Error Type")
ggsave(snakemake@output[["class_by_count_short"]], width = 7, height = 5)

df_vcf_rid_anno_short %>%
  ggplot(aes(y = fct_rev(variant_class), fill = ng_group)) +
  geom_bar() +
  facet_wrap(c("group_name")) +
  labs(x = NULL, y = NULL, fill = "Error count\n(pat, mat)")
ggsave(snakemake@output[["class_by_number_short"]], width = 7, height = 5)

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
  labs(x = "count", y = NULL, fill = "category")
ggsave(snakemake@output[["ave_collapse"]])
