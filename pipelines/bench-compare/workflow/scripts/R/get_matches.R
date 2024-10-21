library(tidyverse)

# Remove variants for which genome is not expected. We can do this iif for any
# variants for which a genome is not expected, all variants within the same
# group for which a genome is expected have a provable match to a genome error
# (either 1-1 or in combination with multiple variants/errors). If this is not
# true, we should not remove a variant that does not have an expected hit since
# this might be necessary to replay the variant to an equivalent error, and the
# reason we cannot prove this is because the genome error lifted over to the
# wrong place.

# In addition, we don't need to keep variants which were duplicated in the
# intersection that have one and only one rowmatch. In this case, we know which
# error matches the variant, so we don't need to keep the other rows in the
# intersection b/t the wrong variant/error.

ensure <- function(test, msg) {
  if (!test) {
    print(test)
    stop(sprintf("ERROR %s", msg))
  }
}

ensure_empty <- function(df, msg) {
  ensure(nrow(df) == 0, msg)
}

bool_comb <- function(n) {
  if (n == 1) {
    matrix(c(FALSE, TRUE), ncol = 1)
  } else {
    w <- 2 ^ (n - 1)
    matrix(
      c(
        rep(FALSE, w),
        rep(TRUE, w),
        matrix(rep(t(bool_comb(n - 1)), 2), ncol = n - 1, byrow = TRUE)
      ),
      ncol = n
    )
  }
}

# this doesn't exist, womp womp :(
group_indices_ordered <- function(grouped_df) {
  grouped_df %>%
    group_indices() %>%
    match(., unique(.))
}

index_group <- function(df, c) {
  add_column(df, {{ c }} := group_indices_ordered(df))
}

read_stuff <- function(path) {
  read_tsv(
    path,
    col_types = cols(
      regions = "c",
      nocov = "-",
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
      trim_left = "-",
      trim_right = "-",
      error_group = "-",
      error_group_size = "-",
      truth_len = "i",
      query_len = "i",
      .default = "c",
      ),
    na = c(".", "NA")
  )
}

replay_variant <- function(refseq, start, end, variant) {
  str_sub(refseq, pmax(0, start), end) <- variant
  refseq
}

get_rowmatch <- function(seq, truth_seq, query_seq, q100_seq, hprc_seq,
                         seqstart, g0, g1, v0, v1) {
  .g0 <- g0 - seqstart + 1
  .g1 <- g1 - seqstart
  .v0 <- v0 - seqstart + 1
  .v1 <- v1 - seqstart
  .q100_rpy <- replay_variant(seq, .g0, .g1, q100_seq)
  .hprc_rpy <- replay_variant(seq, .g0, .g1, hprc_seq)
  .truth_rpy <- replay_variant(seq, .v0, .v1, truth_seq)
  .query_rpy <- replay_variant(seq, .v0, .v1, query_seq)
  if_else(is.na(truth_seq), .q100_rpy == seq, .q100_rpy == .truth_rpy) &
    if_else(is.na(query_seq), .hprc_rpy == seq, .hprc_rpy == .query_rpy)
}

#
# read in lots of stuff
#

df_g0 <- read_stuff(snakemake@input[["gbench"]])

df_v0 <- read_stuff(snakemake@input[["vbench"]])

df_v1 <- df_v0 %>%
  filter(!vchrom %in% c("chrX", "chrY")) %>%
  filter(src_hap == dst_hap) %>%
  select(-dst_hap) %>%
  rename(hap = src_hap) %>%
  mutate(
    liftover_diff = q100_len - (realend - realstart),
    across(
      c(ref, truth_alt, query_alt, q100, hprc),
      ~ if_else(.x == "*", "", .x)
    ),
    chrom = as.integer(str_sub(vchrom, 4)),
  ) %>%
  arrange(desc(hap), vid, gid) %>%
  select(-alt, -q100_len, -hprc_len, -truth_len, -query_len, -vchrom, -gchrom) %>%
  relocate(regions, .after = last_col())

#
# Remove genome errors with obvious liftover errors (difference >=100 bp)
#

df_v_badlift <- df_v1 %>%
  filter(!(abs(liftover_diff) < 100 | is.na(gstart)))

df_v <- df_v1 %>%
  anti_join(df_v_badlift, by = c("hap", "vid")) %>%
  relocate(vid, gid) %>%
  select(-liftover_diff)

#
# Clean up variant/error hits and find 1-1 matches
#
# Steps to do this:
# 1. remove exact duplicate errors (we can only use one of them)
# 2. attempt to find 1-1 matches
# 3. fix overlapping genome errors (easier to do with 1-1 matches from (2))
# 4. attempt to find more 1-1 matches after fixing overlaps
# 5. assign all rows into "intersection groups"

# These are genome errors that exactly overlap another genome error and encode
# the same edit. It only makes sense to keep one of these, so ignore the rest
# but keep them around so we can examine them if we want.
df_v_g_duplicated <- df_v %>%
  filter(!is.na(gid)) %>%
  group_by(hap, chrom, realstart, realend, q100, hprc) %>%
  filter(gid != first(gid)) %>%
  ungroup() %>%
  select(hap, gid, chrom, realstart, realend) %>%
  unique()

df_v_rowmatch <- df_v %>%
  filter(!is.na(gid)) %>%
  # Remove the duplicates we found above
  anti_join(df_v_g_duplicated, by = c("hap", "gid")) %>%
  mutate(
    # Compute the length differences of each edit; use this later to match cases
    # where the liftover "almost worked". 
    ref_len = str_length(ref),
    truth_indel_len = replace_na(str_length(truth_alt) - ref_len, 0),
    query_indel_len = replace_na(str_length(query_alt) - ref_len, 0),
    .lift_len = realend - realstart,
    .q100_indel_len = str_length(q100) - .lift_len,
    .hprc_indel_len = str_length(hprc) - .lift_len,
    .truth_diff = truth_indel_len - .q100_indel_len,
    .query_diff = query_indel_len - .hprc_indel_len,
    # Determine which rows "almost" match 1-1 and will match if we move the
    # starting position of the liftover (because it failed to land in the
    # proper spot); diff equality check is to ensure that we only match if
    # the edit lengths between truth and query are the same.
    single_adj_realstart = realstart + .truth_diff,
    rowmatch = .truth_diff == .query_diff &
      get_rowmatch(seq, truth_alt, query_alt, q100, hprc, gstart,
                   single_adj_realstart, realend, vstart, vend)
  ) %>%
  # ASSUME at this point, if a vid matches a gid 1-1, then the gid/vid match
  # *nowhere else* (ie there are no duplicates...because we removed them).
  # Standardize the boundaries for each gid so they reflect the adjustments made
  # above.
  group_by(hap, gid) %>%
  mutate(
    .s = first(single_adj_realstart[rowmatch]),
    # revert back to original realstart for non-matches since we have no
    # evidence upon which to base our adjustment
    single_adj_realstart = if_else(is.na(.s), realstart, .s)
  ) %>%
  ungroup() %>%
  select(-starts_with("."))

# Some genome errors overlap each other exactly but don't encode the same edit.
# Many of these cases are pairs of SNP and INDEL genome errors that are adjacent
# to each other. In these cases, it is usually enough to move the start of the
# SNP error such that the error is 1bp (ie the end doesn't move) and then adjust
# the end of the INDEL error to the start of the SNP. The beginning of the INDEL
# may need to be further adjusted using the rowmatch logic described above.
#
# There are a few overlaps that are not like this that we can ignore for now
# since they are so few.
df_v_perfect_overlaps <- df_v_rowmatch %>%
  mutate(.shift = single_adj_realstart - realstart) %>%
  group_by(gid, hap, chrom, realstart, realend, .shift, q100, hprc) %>%
  summarize(.any_rowmatch = any(rowmatch), .groups = "drop") %>%
  group_by(hap, chrom, realstart, realend) %>%
  index_group(overlap_id) %>%
  group_by(overlap_id) %>%
  mutate(
    overlap_n = n(),
    .is_snp = str_length(q100) == 1 & str_length(hprc) == 1
  ) %>%
  mutate(is_snp_overlap = overlap_n == 2 & sum(.is_snp) >= 1) %>%
  ungroup() %>%
  filter(overlap_n > 1) %>%
  arrange(overlap_id, .is_snp, .any_rowmatch) %>%
  group_by(overlap_id) %>%
  mutate(.is_last_gid = last(gid) == gid) %>%
  ungroup() %>%
  mutate(
    overlap_realstart =
      case_when(
        is_snp_overlap & .is_last_gid ~ realstart + .shift,
        is_snp_overlap ~ realstart,
        TRUE ~ NA
      ),
    overlap_realend =
      case_when(
        is_snp_overlap & .is_last_gid ~ realend,
        is_snp_overlap ~ realstart + lead(.shift),
        TRUE ~ NA
      )
  ) %>%
  select(-starts_with("."), -q100, -hprc, -realstart, -realend, -chrom)

df_v_fixed_overlaps <- df_v_perfect_overlaps %>%
  right_join(df_v_rowmatch, by = c("gid", "hap")) %>%
  replace_na(list(is_snp_overlap = FALSE)) %>%
  mutate(
    .lift_len = overlap_realend - overlap_realstart,
    .q100_indel_len = str_length(q100) - .lift_len,
    .hprc_indel_len = str_length(hprc) - .lift_len,
    .truth_diff = truth_indel_len - .q100_indel_len,
    .query_diff = query_indel_len - .hprc_indel_len,
    overlap_adj_realstart = overlap_realstart + .truth_diff,
    overlap_rowmatch = if_else(
      is.na(overlap_id),
      rowmatch,
      get_rowmatch(seq, truth_alt, query_alt, q100, hprc, gstart,
                   overlap_adj_realstart, overlap_realend, vstart, vend)
    )
  ) %>%
  select(-starts_with("."))

# TEST: anything that was a rowmatch is still a rowmatch after fixing the
# overlaps
df_v_fixed_overlaps %>%
  filter(rowmatch & !overlap_rowmatch) %>%
  ensure_empty("All previous rowmatches must still be overlapping rowmatches")

# Combine all the adjustments and matches we just did (or didn't do)
df_v_fixed_rowmatch <- df_v_fixed_overlaps %>%
  mutate(rowmatch = rowmatch | replace_na(overlap_rowmatch, FALSE)) %>%
  select(-overlap_rowmatch) %>%
  mutate(
    # Only keep the adjusted start position if we have evidence that it needed
    # to be moved in the form of a rowmatch (either before or after fixing
    # overlaps). Since we know that each gid has 0 or 1 rowmatches, this also
    # implies that each gid will have 0 or 1 adjusted realstart values
    adj_realstart = case_when(
      rowmatch & is_snp_overlap ~ overlap_adj_realstart,
      rowmatch ~ single_adj_realstart,
      TRUE ~ NA
    ),
    adj_realend = if_else(is_snp_overlap, overlap_realend, realend)
  ) %>%
  group_by(hap, gid) %>%
  # Set all rows with the same gid to the adj realstart value if there is one
  mutate(adj_realstart = first(adj_realstart, na_rm = TRUE)) %>%
  ungroup() %>%
  select(-overlap_adj_realstart, -single_adj_realstart, -overlap_realstart,
         -overlap_realend)

# TEST: all vids match exactly 1 or 0 times
df_v_fixed_rowmatch %>%
  group_by(vid, hap) %>%
  filter(sum(rowmatch) > 1) %>%
  ensure_empty("All vids must have 0 or 1 rowmatches")

# TEST: all gids match exactly 1 or 0 times
df_v_fixed_rowmatch %>%
  group_by(gid, hap) %>%
  filter(sum(rowmatch) > 1) %>%
  ensure_empty("All gids must have 0 or 1 rowmatches")

df_v_isec_pre <- df_v_fixed_rowmatch %>%
  arrange(desc(hap), chrom, vstart) %>%
  mutate(
    .start = pmin(vstart, gstart),
    .end = pmax(vend, gend)
  ) %>%
  group_by(hap, chrom) %>%
  mutate(.isec = cumsum(.start >= lag(cummax(.end), default = TRUE))) %>%
  ungroup() %>%
  group_by(.isec, hap, chrom) %>%
  mutate(
    isec_n = n(),
    isec_nv = length(unique(vid)),
    isec_ng = length(unique(gid))
  ) %>%
  ungroup() %>%
  group_by(hap, chrom, .isec) %>%
  index_group(isec) %>%
  ungroup() %>%
  select(-.isec) %>%
  select(-starts_with(".")) %>%
  relocate(isec, isec_n, isec_nv, isec_ng, hap, gid, vid)

df_v_unfixable_overlaps <- df_v_isec_pre %>%
  group_by(isec) %>%
  filter(sum(!is.na(overlap_id) & !is_snp_overlap) > 0) %>%
  ungroup()

df_v_isec <- df_v_isec_pre %>%
  anti_join(df_v_unfixable_overlaps, by = "isec") %>%
  # take these out since at this point the only overlaps are SNP overlaps
  select(-is_snp_overlap, -overlap_n) %>%
  relocate(overlap_id, .before = regions)

#
# Get single hits that are exact matches
#

df_v_1to1 <- df_v_isec %>%
  filter(isec_n == 1)

# These are genome errors that intersected with a variant for which an error
# was not expected, in which case we can ignore the intersection but keep the
# genome error.
df_v_1to1_nohit_g <- df_v_1to1 %>%
  filter(rowmatch & !genome_expected)

# NOTE this will include variants for which a genome hit is expected and not
# expected. In reality, the latter is vanishingly small. These are likely due to
# imperfections in my definition of an expected hit, and regardless this is
# likely not a problem due to the small number
df_v_1to1_matched <- df_v_1to1 %>%
  anti_join(df_v_1to1_nohit_g, by = "isec") %>%
  mutate(
    match_group = if_else(rowmatch, 1, NA),
    match_ng = if_else(rowmatch, 1, NA),
    match_nv = if_else(rowmatch, 1, NA)
  )

#
# clean up multiple variants
#

df_v_many <- df_v_isec %>%
  filter(isec_n > 1)

df_v_many_actions <- df_v_many %>%
  group_by(isec, vid) %>%
  mutate(.v_has_rowmatch = sum(rowmatch) > 0) %>%
  group_by(isec, gid) %>%
  mutate(.g_has_rowmatch = sum(rowmatch) > 0) %>%
  group_by(isec, vid) %>%
  mutate(
    .v_can_kill = !rowmatch & !.v_has_rowmatch & .g_has_rowmatch,
    .v_is_nohit = n() == sum(.v_can_kill)
  ) %>%
  group_by(isec, gid) %>%
  mutate(
    .g_can_kill = !rowmatch & .v_has_rowmatch & !.g_has_rowmatch,
    .g_is_nohit = n() == sum(.g_can_kill)
  ) %>%
  ungroup() %>%
  mutate(
    .action = case_when(
      rowmatch & .v_has_rowmatch & .g_has_rowmatch ~ "single",
      !rowmatch & .v_has_rowmatch & .g_has_rowmatch ~ "kill",
      .g_can_kill ~ if_else(.g_is_nohit, "g_nohit", "kill"),
      .v_can_kill ~ if_else(.v_is_nohit & genome_expected, "v_nohit", "kill"),
      TRUE ~ "combi"
    )
  ) %>%
  group_by(isec) %>%
  mutate(
    action = if_else(
      .action != "combi",
      .action,
      case_when(
        sum(.action == "combi") == 1 ~ if_else(genome_expected, "mismatch", "g_nohit"),
        TRUE ~ .action
      )
    )
  ) %>%
  ungroup() %>%
  select(-starts_with("."))

df_v_many_combi <- df_v_many_actions %>%
  group_by(isec) %>%
  filter(sum(action == "combi") > 0) %>%
  ungroup()

df_v_many_g_nohit <- df_v_many_actions %>%
  anti_join(df_v_many_combi, by = "isec") %>%
  filter(action == "g_nohit")

df_v_many_v_nohit <- df_v_many_actions %>%
  anti_join(df_v_many_combi, by = "isec") %>%
  filter(action == "v_nohit")

df_v_many_matched <- df_v_many_actions %>%
  anti_join(df_v_many_combi, by = "isec") %>%
  anti_join(df_v_many_v_nohit, by = "isec") %>%
  anti_join(df_v_many_g_nohit, by = "isec") %>%
  filter(action != "kill") %>%
  mutate(
    match_group = if_else(action == "single", 1, NA),
    match_ng = if_else(action == "single", 1, NA),
    match_nv = if_else(action == "single", 1, NA)
  )

#
# get multiple variant hits (one g to many v)
#

df_v_1toN <- df_v_many_combi %>%
  filter(isec_ng == 1)

# TODO this currently does not adjust the genome error start position given the
# variant we are trying to match. This would be easy to implement and should
# give (slightly) higher yield
replay_multi <- function(df, key) {
  nvar <- nrow(df)
  seq <- df$seq[[1]]
  q100_seq <- df$.q100_rpy[[1]]
  hprc_seq <- df$.hprc_rpy[[1]]
  v0 <- df$.v0
  v1 <- df$.v1
  truth_alt <- df$truth_alt
  query_alt <- df$query_alt
  cs <- bool_comb(nvar)
  ncomb <- nrow(cs)
  matches <- rep(FALSE, nvar)
  # skip first combination which is entirely FALSE (which will replay nothing)
  for (i in 2:ncomb) {
    truth_seq <- seq
    query_seq <- seq
    # loop through each variant bench line in reverse so that the position
    # doesn't get screwed up
    for (j in nvar:1) {
      if (cs[i, j]) {
        if (!is.na(truth_alt[[j]])) {
          str_sub(truth_seq, v0[[j]], v1[[j]]) <- truth_alt[[j]]
        }
        if (!is.na(query_alt[[j]])) {
          str_sub(query_seq, v0[[j]], v1[[j]]) <- query_alt[[j]]
        }
      }
    }
    if (truth_seq == q100_seq & query_seq == hprc_seq) {
      matches <- cs[i, ]
      break
    }
  }
  matches
  df %>%
    mutate(
      vmatches = matches
    ) %>%
    add_column(
      gid = key$gid[[1]],
      hap = key$hap[[1]]
    )
}

df_v_1toN_results <- df_v_1toN %>%
  group_by(isec) %>%
  mutate(
    # Figure out how much to adjust the liftover region based on edit lengths:
    # this is similar to what we do for single variants (see above) except that
    # a) the sum of all edit length in truth and query should equal the edit
    # length of the genome error (since there is only one) and b) the only
    # variants that count toward this sum are variants that are not TP/TP that
    # don't expect a genome hit
    .is_real_tp = replace_na(truth_bd == "TP" & query_bd == "TP" & !genome_expected, FALSE),
    .lift_len = realend - realstart,
    .q100_ind_len = str_length(q100) - .lift_len,
    .hprc_ind_len = str_length(hprc) - .lift_len,
    .truth_diff = sum(if_else(.is_real_tp, 0, truth_indel_len)) - .q100_ind_len,
    .query_diff = sum(if_else(.is_real_tp, 0, query_indel_len)) - .hprc_ind_len,
    .adj_realstart = if_else(
      .truth_diff == .query_diff & .truth_diff != 0,
      realstart + .truth_diff,
      realstart
    )
  ) %>%
  ungroup() %>%
  # precompute lots of other stuff
  mutate(
    .v0 = pmax(vstart - gstart + 1, 0),
    .v1 = vend - gstart,
    .q100_rpy = replay_variant(seq, .adj_realstart - gstart + 1, realend - gstart, q100),
    .hprc_rpy = replay_variant(seq, .adj_realstart - gstart + 1, realend - gstart, hprc)
  ) %>%
  # ensure variant bench lines are in order
  arrange(hap, gid, chrom, vstart, vend) %>%
  group_by(gid, hap) %>%
  group_map(replay_multi) %>%
  bind_rows() %>%
  select(-starts_with("."))

df_v_1toN_matched <- df_v_1toN_results %>%
  group_by(isec) %>%
  mutate(
    .has_match = sum(vmatches) > 0,
    # count rows where a genome hit was expected but no match found
    .n_missing_matches = sum(genome_expected & !vmatches)
  ) %>%
  # remove rows that have no expected hit and no match but only if they are part
  # of a group that isn't missing a match
  filter(
    .n_missing_matches > 0 |
      (.n_missing_matches == 0 & !(!genome_expected & !vmatches))
  ) %>%
  mutate(
    match_group = if_else(vmatches, 1, NA),
    match_ng = if_else(vmatches, 1, NA),
    match_nv = if_else(vmatches, sum(vmatches), NA)
  ) %>%
  ungroup() %>%
  arrange(isec) %>%
  select(-vmatches)

#
# get multi-hits (many g to many v)
#

# Keep 1-1 rowmatches in this dataset since these might be required to match a
# neighboring error correctly. For example, if an imperfect repeat on hg38 is
# flanked by two variants which match one error at the end of the repeat and
# result in an expansion or contraction, and the middle of the repeat has a 1-1
# variant/error match which makes the repeat perfect. The former match is
# dependent on the latter match to be replayed since otherwise the imperfection
# in the middle of the repeat could end up in a different relative position.

df_v_NtoN <- df_v_many_combi %>%
  filter(isec_ng > 1) %>%
  # For groups that have a single genome error without a direct 1-1 match, the
  # remaining non-matches can be thought as a a sub one-to-many group, and thus
  # the starting position of the one genome error can be fixed like we did above
  group_by(isec) %>%
  mutate(.missing_one_g = sum(rowmatch) == isec_ng - 1) %>%
  group_by(gid) %>%
  mutate(
    .adj_this_g = if_else(.missing_one_g, sum(rowmatch) == 1, FALSE),
    .is_real_tp = replace_na(truth_bd == "TP" & query_bd == "TP" & !genome_expected, FALSE),
    .lift_len = realend - realstart,
    .q100_ind_len = str_length(q100) - .lift_len,
    .hprc_ind_len = str_length(hprc) - .lift_len,
    .truth_diff = sum(if_else(.is_real_tp, 0, truth_indel_len)) - .q100_ind_len,
    .query_diff = sum(if_else(.is_real_tp, 0, query_indel_len)) - .hprc_ind_len,
    adj_realstart = if_else(
      .truth_diff == .query_diff & .truth_diff != 0 & .adj_this_g,
      realstart + .truth_diff,
      adj_realstart
    )
  ) %>%
  ungroup() %>%
  select(-starts_with("."))

df_v_NtoN_notworthit <- df_v_NtoN %>%
  filter((2 ^ isec_nv - 1) * (2 ^ isec_ng - 1) > 1e5)

df_v_NtoN_worthit <- df_v_NtoN %>%
  anti_join(df_v_NtoN_notworthit, by = "isec")

vallele <- setClass(
  "vallele",
  slots = list(
    id = "numeric",
    start = "numeric",
    end = "numeric",
    truth = "character",
    query = "character",
    rowmatch = "logical"
  )
)

gallele <- setClass(
  "gallele",
  slots = list(
    id = "numeric",
    start = "numeric",
    end = "numeric",
    q100 = "character",
    hprc = "character",
    rowmatch = "logical"
  )
)

df_v_poly_seqs <- df_v_NtoN_worthit %>%
  select(isec, seq, gstart) %>%
  unique() %>%
  arrange(isec, gstart) %>%
  rename(seqstart = gstart) %>%
  group_by(isec) %>%
  mutate(
    diff = lead(seqstart) - seqstart,
    seq = if_else(is.na(diff), seq, str_sub(seq, 0, diff))
  ) %>%
  summarize(seq = str_flatten(seq), seqstart = min(seqstart), .groups = "drop")

df_v_valleles <- df_v_NtoN_worthit %>%
  group_by(vid) %>%
  mutate(.has_rowmatch = sum(rowmatch) > 0) %>%
  ungroup() %>%
  select(isec, vid, vstart, vend, truth_alt, query_alt, .has_rowmatch) %>%
  unique() %>%
  left_join(df_v_poly_seqs, by = "isec") %>%
  mutate(
    v0 = pmax(vstart - seqstart + 1, 0),
    v1 = vend - seqstart,
    vallele = pmap(
      list(vid, v0, v1, truth_alt, query_alt, .has_rowmatch),
      function(i, s, e, t, q, r) {
        vallele(id = i, start = s, end = e, truth = t, query = q, rowmatch = r)
      }
    )
  ) %>%
  group_by(isec) %>%
  summarize(
    valleles = list(vallele),
    .groups = "drop"
  )

df_v_galleles <- df_v_NtoN_worthit %>%
  group_by(gid) %>%
  mutate(
    .has_rowmatch = sum(rowmatch) > 0,
    .start = coalesce(adj_realstart, realstart),
    .end = coalesce(adj_realend, realend)
  ) %>%
  ungroup() %>%
  select(isec, gid, .start, .end, q100, hprc, .has_rowmatch) %>%
  unique() %>%
  left_join(df_v_poly_seqs, by = "isec") %>%
  mutate(
    g0 = pmax(.start - seqstart + 1, 0),
    g1 = .end - seqstart,
    gallele = pmap(
      list(gid, g0, g1, q100, hprc, .has_rowmatch),
      function(i, s, e, t, q, r) {
        gallele(id = i, start = s, end = e, q100 = t, hprc = q, rowmatch = r)
      }
    )
  ) %>%
  group_by(isec) %>%
  summarize(
    galleles = list(gallele),
    .groups = "drop"
  )

group_columns <- function(m, indices_only) {
  c <- ncol(m)
  skip <- rep(FALSE, c)
  acc <- list()
  for (i in 1:c) {
    if (!skip[[i]]) {
      match_runs <- m[, i]
      matches <- apply(m == m[, i], 2, all)
      skip <- skip | matches
      if (indices_only) {
        acc[[length(acc) + 1]] <- which(matches)
      } else {
        acc[[length(acc) + 1]] <- list(indices = which(matches), runs = match_runs)
      }
    }
  }
  acc
}

merge_columns <- function(m) {
  c <- ncol(m)
  acc <- replicate(c, NULL, simplify = FALSE)
  for (i in 1:c) {
    if (i < c) {
      for (j in (i+1):c) {
        s0 <- m[, i]
        s1 <- m[, j]
        s01 <- s0 & s1
        if (all(s01 == s0) || all(s01 == s1)) {
          acc[[i]] <- rbind(acc[[i]], s01)
          acc[[j]] <- rbind(acc[[j]], s01)
        }
      }
    }
  }
  acc <- map(acc, unique)
  # I don't know what to do if this test doesn't pass :/
  ensure(every(acc, ~ is.null(.x) || nrow(.x) == 1), "Non-unique merge detected")
  for (i in 1:c) {
    if (!is.null(acc[[i]])) {
      m[, i] <- acc[[i]]
    }
  }
  group_columns(m, TRUE)
}

replay_multi_poly <- function(i, v, g, seq) {
  nv <- length(v)
  ng <- length(g)
  # TODO optimize the combinations we need to search. We know that some
  # combinations match 1-1, so we can filter these matrices to only use the
  # combinations that include these known 1-1 matches. This assumes there are
  # no "double matches" (which should be removed anyway)
  vcombs <- bool_comb(nv)
  gcombs <- bool_comb(ng)
  nvcomb <- nrow(vcombs)
  ngcomb <- nrow(gcombs)
  vres <- NULL
  gres <- NULL
  # skip first combination in both matrices which is entirely FALSE (which will
  # replay nothing)
  for (vi in 2:nvcomb) {
    for (gi in 2:ngcomb) {
      truth_seq <- seq
      query_seq <- seq
      q100_seq <- seq
      hprc_seq <- seq
      # loop through each variant bench line in reverse so that the position
      # doesn't get screwed up
      for (j in nv:1) {
        if (vcombs[vi, j]) {
          if (!is.na(v[[j]]@truth)) {
            str_sub(truth_seq, v[[j]]@start, v[[j]]@end) <- v[[j]]@truth
          }
          if (!is.na(v[[j]]@query)) {
            str_sub(query_seq, v[[j]]@start, v[[j]]@end) <- v[[j]]@query
          }
        }
      }
      # ditto genome bench line
      for (j in ng:1) {
        if (gcombs[gi, j]) {
          if (!is.na(g[[j]]@q100)) {
            str_sub(q100_seq, g[[j]]@start, g[[j]]@end) <- g[[j]]@q100
          }
          if (!is.na(g[[j]]@hprc)) {
            str_sub(hprc_seq, g[[j]]@start, g[[j]]@end) <- g[[j]]@hprc
          }
        }
      }
      ## print(i)
      ## print(truth_seq)
      ## print(q100_seq)
      ## print(query_seq)
      ## print(hprc_seq)
      if (truth_seq == q100_seq & query_seq == hprc_seq) {
        vres <- rbind(vres, vcombs[vi, ])
        gres <- rbind(gres, gcombs[gi, ])
      }
    }
  }
  if (!is.null(vres)) {
    vids <- map(v, ~ .x@id)
    gids <- map(g, ~ .x@id)
    matches <- group_columns(cbind(vres, gres), FALSE)
    match_indices <- map(matches, ~ .x$indices)
    matches %>%
      map(~ .x$runs) %>%
      do.call(cbind, .) %>%
      merge_columns %>%
      map(~ unlist(match_indices[.x])) %>%
      map(~ list(vid = vids[.x[.x <= nv]], gid = gids[.x[.x > nv] - nv]))
  }
}

df_v_NtoN_results <- df_v_poly_seqs %>%
  left_join(df_v_valleles, by = "isec") %>%
  left_join(df_v_galleles, by = "isec") %>%
  select(-seqstart) %>%
  ## head(n = 20) %>%
  ## filter(isec %in% c(239, 1136)) %>%
  ## filter(isec %in% c(32)) %>%
  mutate(result = pmap(list(isec, valleles, galleles, seq), replay_multi_poly)) %>%
  select(-seq, -valleles, -galleles)

df_v_NtoN_matched <- df_v_NtoN_results %>%
  unnest_longer(result) %>%
  unnest_wider(result) %>%
  mutate(
    match_group = row_number(),
    match_ng = map_int(gid, length),
    match_nv = map_int(vid, length)
  ) %>%
  mutate(c = map2(vid, gid, ~ expand_grid(vid = unlist(.x), gid = unlist(.y)))) %>%
  select(-vid, -gid) %>%
  unnest(c) %>%
  right_join(df_v_NtoN_worthit, by = c("isec", "gid", "vid")) %>%
  arrange(isec) %>%
  # if a gid/vid has a match, remove all rows that don't have a match
  group_by(hap, gid) %>%
  mutate(.g_has_match = sum(!is.na(match_group)) > 0) %>%
  filter(if_else(.g_has_match, !is.na(match_group), is.na(match_group))) %>%
  group_by(hap, vid) %>%
  mutate(.v_has_match = sum(!is.na(match_group)) > 0) %>%
  filter(if_else(.v_has_match, !is.na(match_group), is.na(match_group))) %>%
  ungroup() %>%
  select(-starts_with(".")) %>%
  bind_rows(df_v_NtoN_notworthit)

#
# output everything
#

# all match/unmatch results
df_v_all_matched <- bind_rows(
  df_v_1to1_matched,
  df_v_many_matched,
  df_v_1toN_matched,
  df_v_NtoN_matched
) %>%
  arrange(isec, vstart) %>%
  select(-starts_with("."), -ends_with("_len"), -gstart, -gend, -seq) %>%
  group_by(isec, match_group) %>%
  index_group(match_id) %>%
  ungroup() %>%
  mutate(match_id = if_else(!is.na(match_group), match_id, NA)) %>%
  select(-match_group) %>%
  relocate(regions, .after = last_col()) %>%
  relocate(match_id, .before = match_ng) %>%
  write_tsv(snakemake@output[["hits"]])
  
## df_v_unmatched_nohits %>%
##   mutate(vchrom = sprintf("chr%d", vchrom)) %>%
##   mutate(name = sprintf("%s->%s", truth_alt, query_alt)) %>%
##   relocate(vchrom, vstart, vend, name) %>%
##   select(vchrom, vstart, vend, name, hap, vid) %>%
##   write_tsv("vbench_nohits.bed", col_names = F)

df_v %>%
  filter(is.na(gid) & genome_expected) %>%
  bind_rows(df_v_many_v_nohit) %>%
  select(-gstart, -gend, -realstart, -realend, -gid, -error_type,
         -q100, -hprc, -seq, -ends_with("_len")) %>%
  relocate(chrom, vstart, vend) %>%
  write_tsv(snakemake@output[["v_nohit"]])

df_g0 %>%
  filter(is.na(vid)) %>%
  filter(src_hap == dst_hap) %>%
  rename(hap = src_hap) %>%
  select(-dst_hap) %>%
  filter(!gchrom %in% c("chrX", "chrY")) %>%
  mutate(chrom = as.integer(str_sub(gchrom, 4))) %>%
  bind_rows(df_v_1to1_nohit_g) %>%
  bind_rows(df_v_many_g_nohit) %>%
  arrange(gid, desc(hap)) %>%
  select(-vchrom, -gchrom, -vstart, -vend, -vid, -ref, -alt, -regions, -starts_with("truth"),
         -starts_with("query"), -genome_expected, -seq, -ends_with("_len")) %>%
  relocate(chrom) %>%
  write_tsv(snakemake@output[["g_nohit"]])

df_v_badlift %>%
  write_tsv(snakemake@output[["badlift"]])

df_v_unfixable_overlaps %>%
  write_tsv(snakemake@output[["unfixable_overlap"]])

df_v_g_duplicated %>%
  write_tsv(snakemake@output[["duplicated"]])
