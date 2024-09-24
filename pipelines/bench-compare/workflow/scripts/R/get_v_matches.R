library(tidyverse)

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

df_v_pre <- read_stuff(snakemake@input[[1]]) %>%
  filter(!vchrom %in% c("chrX", "chrY")) %>%
  mutate(vchrom = as.integer(str_sub(vchrom, 4))) %>%
  arrange(vid, desc(hap)) %>%
  mutate(
    hasg = !is.na(gchrom),
    ng = if_else(hasg, ng, 0)
  )

df_v_badlift <- df_v_pre %>%
  filter(!(abs(liftover_diff) < 100 | is.na(gchrom)))
  
df_v_badlift %>%
  write_tsv(snakemake@output[["badlift"]])

df_v <- df_v_pre %>%
  anti_join(df_v_badlift, by = c("hap", "vid"))
 
#
# get single hits that are exact matches
#

replay_variant <- function(refseq, start, end, variant) {
  str_sub(refseq, pmax(0, start), end) <- variant
  refseq
}

df_v_1_replay <- df_v %>%
  filter(hasg) %>%
  filter(nv == 1) %>%
  filter(ng == 1) %>%
  mutate(
    q100_rpy = replay_variant(seq, realstart - gstart + 1, realend - gstart, q100),
    hprc_rpy = replay_variant(seq, realstart - gstart + 1, realend - gstart, hprc),
    truth_rpy = replay_variant(seq, vstart - gstart + 1, vend - gstart, truth_alt),
    query_rpy = replay_variant(seq, vstart - gstart + 1, vend - gstart, query_alt),
  )

df_v_1_matched <- df_v_1_replay %>%
  filter(
    if_else(is.na(truth_alt), q100_rpy == seq, q100_rpy == truth_rpy)
    & if_else(is.na(query_alt), hprc_rpy == seq, hprc_rpy == query_rpy))

df_v_1_unmatched <- df_v_1_replay %>%
  anti_join(df_v_1_matched, by = c("vid", "hap"))

#
# get multiple variant hits (one g to many v)
#

df_v_ggroup <- df_v %>%
  filter(ng > 1) %>%
  filter(nv == 1) %>%
  arrange(gid, hap, vid) %>%
  relocate(gid, vid, hap, ng) %>%
  group_by(hap, gid) %>%
  filter(n() == ng) %>%
  ungroup()

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
  df %>%
    mutate(
      vmatches = matches
    ) %>%
    add_column(
      gid = key$gid[[1]],
      hap = key$hap[[1]]
    )
}

df_v_multi <- df_v_ggroup %>%
  # precompute lots of stuff
  mutate(
    .v0 = pmax(vstart - gstart + 1, 0),
    .v1 = vend - gstart,
    .q100_rpy = replay_variant(seq, realstart - gstart + 1, realend - gstart, q100),
    .hprc_rpy = replay_variant(seq, realstart - gstart + 1, realend - gstart, hprc)
  ) %>%
  # ensure variant bench lines are in order
  arrange(hap, gid, vchrom, vstart, vend) %>%
  group_by(gid, hap) %>%
  group_map(replay_multi) %>%
  bind_rows() %>%
  select(-starts_with("."))

# All variant bench lines here have been shown to either be part of a genome
# bench variant or not. Those outside this set (for which no matches were found
# at all) are "inconclusive, and likely failed due to liftover or something.
df_v_multi_matched <- df_v_multi %>%
  group_by(gid, hap) %>%
  filter(sum(vmatches) > 0) %>%
  ungroup()

#
# get multi-hits (many g to many v)
#

df_v_ggroup_poly <- df_v %>%
  filter(ng >= 1) %>%
  filter(nv > 1) %>%
  arrange(gid, vid, hap) %>%
  relocate(gid, vid, hap, ng, nv) %>%
  group_by(hap, gid) %>%
  filter(n() == ng) %>%
  group_by(hap, vid) %>%
  filter(n() == nv) %>%
  ungroup() %>%
  group_by(hap) %>%
  mutate(
    .group = cumsum(
      gid > lag(cummax(gid), default = TRUE) &
        vid > lag(cummax(vid), default = TRUE)
     )
  ) %>%
  ungroup()

vallele <- setClass(
  "vallele",
  slots = list(
    id = "numeric",
    start = "numeric",
    end = "numeric",
    truth = "character",
    query = "character"
  )
)

gallele <- setClass(
  "gallele",
  slots = list(
    id = "numeric",
    start = "numeric",
    end = "numeric",
    q100 = "character",
    hprc = "character"
  )
)

df_v_poly_seqs <- df_v_ggroup_poly %>%
  select(hap, .group, seq, gstart) %>%
  unique() %>%
  arrange(hap, .group) %>%
  rename(seqstart = gstart) %>%
  group_by(hap, .group) %>%
  mutate(
    diff = lead(seqstart) - seqstart,
    seq = if_else(is.na(diff), seq, str_sub(seq, 0, diff))
  ) %>%
  summarize(seq = str_flatten(seq), seqstart = min(seqstart), .groups = "drop")
  
df_v_valleles <- df_v_ggroup_poly %>%
  select(hap, .group, vid, vstart, vend, truth_alt, query_alt) %>%
  unique() %>%
  left_join(df_v_poly_seqs, by = c("hap", ".group")) %>%
  mutate(
    v0 = pmax(vstart - seqstart + 1, 0),
    v1 = vend - seqstart,
    vallele = pmap(
      list(vid, v0, v1, truth_alt, query_alt),
      function(i, s, e, t, q) {
        vallele(id = i, start = s, end = e, truth = t, query = q)
      }
    )
  ) %>%
  group_by(hap, .group) %>%
  summarize(
    valleles = list(vallele),
    .groups = "drop"
  )

df_v_galleles <- df_v_ggroup_poly %>%
  select(hap, .group, gid, realstart, realend, q100, hprc) %>%
  unique() %>%
  left_join(df_v_poly_seqs, by = c("hap", ".group")) %>%
  mutate(
    g0 = pmax(realstart - seqstart + 1, 0),
    g1 = realend - seqstart,
    gallele = pmap(
      list(gid, g0, g1, q100, hprc),
      function(i, s, e, t, q) {
        gallele(id = i, start = s, end = e, q100 = t, hprc = q)
      }
    )
  ) %>%
  group_by(hap, .group) %>%
  summarize(
    galleles = list(gallele),
    .groups = "drop"
  )

replay_multi_poly <- function(v, g, seq) {
  nv <- length(v)
  ng <- length(g)
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
      if (truth_seq == q100_seq & query_seq == hprc_seq) {
        .vres <- map(v[vcombs[vi, ]], ~ .x@id)
        .gres <- map(g[gcombs[gi, ]], ~ .x@id)
        if (length(.vres) > length(vres) && length(.gres) > length(gres)) {
          vres <- .vres
          gres <- .gres
        }
      }
    }
  }
  list(
    vmatches = vres,
    gmatches = gres
  )
}

df_v_poly_results <- df_v_poly_seqs %>%
  left_join(df_v_valleles, by = c("hap", ".group")) %>%
  left_join(df_v_galleles, by = c("hap", ".group")) %>%
  select(-seqstart) %>%
  mutate(result = pmap(list(valleles, galleles, seq), replay_multi_poly)) %>%
  mutate(
    .nv = map_int(valleles, length),
    .ng = map_int(galleles, length),
    vids = map(result, ~ .x[[1]]),
    gids = map(result, ~ .x[[2]]),
    nv = map_int(vids, length),
    ng = map_int(gids, length)
  ) %>%
  ## select(-result, -valleles, -galleles, -seq)
  select(-result, -seq)

df_v_poly_vkeep <- df_v_poly_results %>%
  select(hap, .group, vids) %>%
  unnest_longer(vids, values_to = "vid") %>%
  add_column(vkeep = TRUE)

df_v_poly_gkeep <- df_v_poly_results %>%
  select(hap, .group, gids) %>%
  unnest_longer(gids, values_to = "gid") %>%
  add_column(gkeep = TRUE)

df_v_poly_matched <- df_v_ggroup_poly %>%
  left_join(df_v_poly_vkeep, by = c("vid", "hap", ".group")) %>%
  left_join(df_v_poly_gkeep, by = c("gid", "hap", ".group")) %>%
  filter(gkeep & vkeep) %>%
  select(-gkeep, -vkeep, -.group)

#
# output stuff everything
#

df_v_matched <- bind_rows(
  mutate(df_v_1_matched, match_type = "singleton"),
  mutate(df_v_poly_matched, match1 = "many-many"),
  mutate(df_v_multi_matched, match1 = "one-many")
) %>%
  arrange(vchrom, vstart, vend)

df_v_matched %>%
  write_tsv(snakemake@output[["matched"]])

df_v_unmatched <- df_v %>%
  anti_join(df_v_matched, by = c("hap", "vid"))

# This is the stuff that has yet to be explained that has at least one genome hit
df_v_unmatched %>%
  filter(ng > 0) %>%
  select(-alt) %>%
  write_tsv(snakemake@output[["unmatched"]])

# These are the things that should have a hit but don't (relatively few)
df_v_unmatched_nohits <- df_v_unmatched %>%
  filter(ng == 0) %>%
  filter(genome_expected) %>%
  select(-alt) %>%
  select(-gchrom, -gstart, -gend, -realstart, -realend, -gid, -error_type,
         -q100, -hprc, -q100_len, -hprc_len, -seq) %>%
  select(-ng) %>%
  select(-liftover_diff, -genome_expected)

df_v_unmatched_nohits %>%
  write_tsv(snakemake@output[["nohit"]])

df_v_unmatched_nohits %>%
  mutate(vchrom = sprintf("chr%d", vchrom)) %>%
  mutate(name = sprintf("%s->%s", truth_alt, query_alt)) %>%
  relocate(vchrom, vstart, vend, name) %>%
  select(vchrom, vstart, vend, name, hap, vid) %>%
  write_tsv(snakemake@output[["nohit_bed"]], col_names = F)
