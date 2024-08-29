library(tidyverse)

ensure <- function(test, msg) {
  if (!test) {
    print(test)
    stop(sprintf("ERROR %s", msg))
  }
}

ensure_empty <- function(df, msg) {
  ensure(nrow(df) == 0, msg)
}

read_data <- function(happy_path, merged_path) {
  # Read both vcf files, add an ID field to the happy vcf to keep track of the
  # original variant since we are gonna slice this dataframe to bits later.
  df_merged <- read_tsv(
    merged_path,
    comment = "#",
    col_types = "ci-cc----ccc",
    col_names = c("chrom", "start", "ref", "alt", "truth", "query", "dipcall")
  )
  df_happy <- read_tsv(
    happy_path,
    comment = "#",
    col_types = "ci-cccccccc",
    col_names = c(
      "chrom_happy",
      "start_happy",
      "ref_happy",
      "alt_happy",
      "qual",
      "filter",
      "info",
      "format",
      "truth_happy",
      "query_happy"
    )
  ) %>%
    mutate(id = row_number())

  # Test: the two vcfs are the same length, if not this will end badly
  ensure(
    nrow(df_merged) == nrow(df_happy),
    "happy and merged vcf must be same length"
  )

  # Combine both vcfs horizontally. Arrange by chrom, start, and truth for each
  # so that variants at the same position have the same order wrt to truth
  # phasing
  df <- bind_cols(
    df_merged %>%
      arrange(chrom, start, truth),
    df_happy %>%
      select(chrom_happy, start_happy, id, ref_happy, alt_happy, truth_happy, query_happy) %>%
      arrange(chrom_happy, start_happy, truth_happy)
  ) %>%
    # repair truth field, for some reason bcftools merge will turn 0|. into 0/.
    mutate(
      truth = if_else(
        str_detect(truth, "^0/\\.") & str_detect(truth_happy, "^0|\\."),
        str_replace(truth, "^0/\\.", "0|."),
        truth
      )
    )

  # Test: all coords should be the same and all truth/query fields should be the same
  df %>%
    filter(!(
      chrom == chrom_happy
      & start == start_happy
      & truth == truth_happy
      & query == query_happy
    )) %>%
    ensure_empty("positions/chromosomes don't line up")

  # Parse out all the genotype fields
  df_gt <- df %>%
    select(-chrom_happy, -start_happy, -truth_happy, -query_happy) %>%
    mutate(
      truth_gt = str_replace(truth, ":.*", ""),
      query_gt = str_replace(query, ":.*", ""),
      dipcall_gt = str_replace(dipcall, ":.*", ""),
      # ASSUME that the only unphased GT in truth are homs so mat/pat make sense
      tpat = as.integer(str_extract(truth_gt, "([0-9])(\\/|\\|)([0-9])", 1)),
      tmat = as.integer(str_extract(truth_gt, "([0-9])(\\/|\\|)([0-9])", 3)),
      q1 = as.integer(str_extract(query_gt, "([0-9])(\\/|\\|)([0-9])", 1)),
      q2 = as.integer(str_extract(query_gt, "([0-9])(\\/|\\|)([0-9])", 3)),
      # ASSUME everything from dipcall is phased so mat/pat make sense
      dpat = if_else(
        dipcall_gt == "1",
        1,
        as.integer(str_extract(dipcall_gt, "([0-9])(\\/|\\|)([0-9])", 1))
      ),
      dmat = if_else(
        dipcall_gt == "1",
        1,
        as.integer(str_extract(dipcall_gt, "([0-9])(\\/|\\|)([0-9])", 3))
      ),
      ) %>%
    select(-dipcall)
  # keep this around so it is easy to add back when done
  df_happy_meta <- df_happy %>%
    select(id, qual, filter, info, format)
  list(
    gt = df_gt,
    meta = df_happy_meta
  )
}

split_final <- function(final, df) {
  list(
    final = final,
    remainder = anti_join(df, final, by = "id")
  )
}

# Fix matching hom and null GTs
fix_matching_homs <- function(df) {
  df %>%
    mutate(.query_gt = case_when(
      query_gt == "1/1" & (dipcall_gt == "1|1" | dipcall_gt == "1") ~ "1|1",
      query_gt == "." & dipcall_gt == "./." ~ ".",
      TRUE ~ NA)
      ) %>%
    filter(!is.na(.query_gt)) %>%
    mutate(query_gt = .query_gt) %>%
    select(-.query_gt) %>%
    split_final(df)
}

# Fix matching het GTs (which assumes the ref and alt fields b/t merged and
# original also match)
fix_matching_hets <- function(df) {
  # if the genotype of the query and dipcall matches, then the ref and alts
  # also match b/t merged and unmerged
  df %>%
    filter((q1 == dpat & q2 == dmat) | (q1 == dmat & q2 == dpat)) %>%
    filter(!(ref == ref_happy & alt == alt_happy)) %>%
    ensure_empty("matching het GTs must also have same ref and alt")
  df %>%
    filter((q1 == dpat & q2 == dmat) | (q1 == dmat & q2 == dpat)) %>%
    mutate(query_gt = dipcall_gt) %>%
    split_final(df)
}

# Remove star alleles completely since these shouldn't be in happy at all
fix_star_alleles <- function(df) {
  df %>%
    filter(str_detect(alt, "\\*$")) %>%
    # ASSUME this only applies to hetalts
    filter(query_gt == "0/1") %>%
    # ASSUME the * allele is always last and that the higher number in the GT
    # field corresponds to this. Thus setting the higher number to 0 will create
    # a hetalt (which should be what we want)
    mutate(
      .query_gt = case_when(
        dipcall_gt == "1|2" ~ "1|0",
        dipcall_gt == "2|1" ~ "0|1",
        dipcall_gt == "1|3" ~ "1|0",
        dipcall_gt == "3|1" ~ "0|1",
        TRUE ~ NA
      ),
      .alt = str_replace(alt, ",\\*$", ""),
      .trimlen = str_length(ref) - str_length(ref_happy),
      .ref = str_sub(ref, 0, -1 - .trimlen),
      .alts = str_split(.alt, ","),
      .alt_trimmed = map2_chr(
        .alts,
        .trimlen,
        function(as, l) str_flatten(map_chr(as, ~ str_sub(.x, 0, -1 - l)), ",")
      )
    ) %>%
    filter(!is.na(.query_gt)) %>%
    mutate(
      query_gt = .query_gt,
      alt = .alt_trimmed,
      ref = .ref
    ) %>%
    select(-starts_with(".")) %>%
    split_final(df)
}

split_maybe <- function(xs) {
  if_else(xs == ".", NA, str_split(xs, ","))
}

remove_maybe <- function(xs, ys, x, y) {
  if (y %in% ys) { xs } else { discard(xs, ~ .x == x) }
}

fix_alt <- function(category, alt_happy, alts, alts_happy, pat, mat,
                    pat_trimmed, mat_trimmed, trimlen) {
  if (category == "missing1" || category == "missing2") {
    alt_happy
  } else {
    x <- remove_maybe(alts, alts_happy, pat, pat_trimmed) %>%
      remove_maybe(alts_happy, mat, mat_trimmed) %>%
      map_chr(~ str_sub(.x, 0, -1 - trimlen))
    if (length(x) == 0) { "." } else { str_flatten(x, ",") }
  }
}

# Represents an extra allele to be moved from one variant to another
extra_allele <- setClass(
  "extra_allele",
  slots = list(
    start = "numeric", # vcf start position
    ref = "character", # vcf ref field
    alt = "character", # the alt to be moved (just one alt, not the entire field)
    ispat = "logical"  # whether the alt is paternal or maternal
  )
)

# Represents a variant that is missing either 1 or 2 alleles
variant_missing <- setClass(
  "variant_missing",
  slots = list(
    start = "numeric",  # vcf start position
    ref = "character",  # vcf ref field
    alts = "character", # vcf alts field (split into a vector)
    qpat = "numeric",   # the paternal query allele index (to be assigned)
    qmat = "numeric",   # the paternal query allele index (to be assigned)
    allow2 = "logical"  # if TRUE this variant requires 2 matches, otherwise 1
  )
)

# helper function to precompute missing variant struct
build_missing_struct <- function(start, ref, alts, allow2) {
  if (is.na(allow2)) {
    NA
  } else {
    variant_missing(
      start = start,
      ref = ref,
      alts = alts,
      qpat = 0,
      qmat = 0,
      allow2 = allow2
    )
  }
}

# helper function to precompute extra allele struct
build_extra_struct <- function(category, start, ref, extra_alt, extra_is_paternal,
                               dpat_alt, dmat_alt, dpat, dmat) {
  # TODO could also code for "extra_hetalt" but these don't seem to exist
  if (category == "extra_het") {
    list(
      extra_allele(start = start, ref = ref, alt = extra_alt, ispat = extra_is_paternal)
    )
  } else if (category == "extra_null") {
    if (dpat > 0 && dmat > 0) {
      list(
        extra_allele(start = start, ref = ref, alt = dpat_alt, ispat = TRUE),
        extra_allele(start = start, ref = ref, alt = dmat_alt, ispat = FALSE)
      )
    } else if (dpat > 0) {
      list(
        extra_allele(start = start, ref = ref, alt = dpat_alt, ispat = TRUE)
        )
    } else if (dmat > 0) {
      list(
        extra_allele(start = start, ref = ref, alt = dmat_alt, ispat = FALSE)
        )
    } else {
      NA
    }
  } else {
    NA
  }
}

match_found <- function(missing) {
  ifelse(
    missing@allow2, 
    missing@qpat > 0 && missing@qmat > 0,
    missing@qpat > 0 || missing@qmat > 0
  )
}

# Loop through the available extra alleles and attempt to match them with a
# variant that needs an allele to phase it. This is done per-group. First try
# to use "exact matching" which will only work if the position/ref/alt fields
# all match exactly. In some cases, a matching allele will be shifted by a few
# bases but will still have the same representation as an extra allele that
# occurred earlier. These require a more sophisticated matching strategy.
fix_missing_alleles <- function(tophase, extra) {
  xlen <- length(extra)
  plen <- length(tophase)
  npat <- 0
  n_unmatched <- 0
  for (x in extra) {
    this_alt <- NA
    match_idx <- 0
    for (i in 1:plen) {
      p <- tophase[[i]]
      if (!match_found(p)) {
        if (x@start == p@start && x@ref == p@ref && x@alt %in% p@alts) {
          # try exact matching first
          match_idx <- i
          this_alt <- x@alt
          break
        } else {
          # else try to shave off parts of the ref/alt depending on starting
          # position and ref length
          pos_diff <- p@start - x@start
          if (pos_diff > 0) {
            ref_diff <- (nchar(x@ref) - pos_diff) - nchar(p@ref)
            .ref <- str_sub(x@ref, 1 + pos_diff, -1 - ref_diff)
            .alt <- str_sub(x@alt, 1 + pos_diff, -1 - ref_diff)
            if (.ref == p@ref && .alt %in% p@alts) {
              this_alt <- .alt
              match_idx <- i
              break
            }
          }
        }
      }
    }
    # ASSUME match_idx is also set if this is NA
    if (!is.na(this_alt)) {
      gt_idx <- which(tophase[[match_idx]]@alts == this_alt)[[1]]
      if (x@ispat) {
        npat <- npat + 1 # for counting pats and mats
        tophase[[match_idx]]@qpat <- gt_idx
      } else {
        tophase[[match_idx]]@qmat <- gt_idx
      }
    } else {
      n_unmatched <- n_unmatched + 1
    }
  }
  nmat <- xlen - npat
  all_phased <- TRUE
  pat <- rep(0, plen)
  mat <- rep(0, plen)
  # Test that all variants with missing variants were assigned an allele. The
  # assumption is that a variant missing one allele will still have two 0s in
  # the pat/mat gt fields, and a variant missing two alleles will have at least
  # one 0 in either field (note that 0 is the starting value and basically means
  # "missing" in this context).
  for (i in 1:plen) {
    p <- tophase[[i]]
    if (!match_found(p)) {
      all_phased <- FALSE
      break
    }
    pat[[i]] <- p@qpat
    mat[[i]] <- p@qmat
  }
  list(
    gt = sprintf("%d|%d", pat, mat),
    # Test that the p/maternal allele counts are balanced (ie we put them
    # in the appropriate field)
    pat_balanced = sum(pat > 0) == npat,
    mat_balanced = sum(mat > 0) == nmat,
    # Test that we don't have any leftover matches
    all_matched = n_unmatched == 0,
    all_phased = all_phased
  )
}

fix_var_groups <- function(df, key) {
  # This is why order matters. We assume that the extra GTs (which are already
  # known) are at the top and we append the missing GTs below in place of the
  # NAs which were once there. This also assumes that the order of the missing
  # structs didn't change in the matching loop when called above. This is purely
  # for speed. If we didn't do this, we would need to filter the input dataframe
  # and recombine it, which is hilariously slow and not at all necessary. It
  # doesn't require strict order though :/
  nextra <- sum(is.na(df$.allow2))
  extra_query_gt <- head(df$.query_gt, n = nextra)
  extra <- flatten(head(df$.extra_struct, n = nextra))
  result <- df$.missing_struct %>%
    tail(n = -nextra) %>%
    fix_missing_alleles(extra)
  mutate(
    df,
    chrom = key$chrom[[1]],
    .query_gt = c(extra_query_gt, result$gt),
    # add lots of "unit test" results so we can ensure this didn't screw up,
    # and make it granular enough that we can debug if we do screw up
    .pat_balanced = result$pat_balanced,
    .mat_balanced = result$mat_balanced,
    .all_matched = result$all_matched,
    .all_phased = result$all_phased
  )
}

pick_alt <- function(ref, alts, x) {
  if (is.na(x)) {
    NA
  } else if (x == 0) {
    ref
  } else {
    alts[[x]]
  }
}

# Special case of fixing grouped variants. All variants here should be balanced
# and will be solved by "swapping" alleles within a group.
fix_balanced_groups <- function(df) {
  # Precompute lots of stuff here so I don't need to clutter the tight looping
  # functions that will iterate through each group
  df_precomputed <- df %>%
    mutate(
      # ASSUME there are no NAs here, will use these to filter each group later
      .category = case_when(
        query_gt == "0/1" & dipcall_gt == "./." ~ "missing1",
        query_gt == "2/1" & dipcall_gt == "./." ~ "missing2",
        query_gt == "." & dipcall_gt != "./." ~ "extra_null",
        query_gt == "0/1" & dipcall_gt != "./." ~ "extra_het",
        query_gt == "2/1" & dipcall_gt != "./." ~ "extra_hetalt"
      ),
      # Used for missing variants only; missing variants are missing either 1
      # or 2 alleles depending on the genotype
      .allow2 = case_when(
        .category == "missing1" ~ FALSE,
        .category == "missing2" ~ TRUE,
        TRUE ~ NA
      ),
      .alts = split_maybe(alt),
      .alts_happy = split_maybe(alt_happy),
      ## .tpat_alt = pmap_chr(list(ref, .alts_happy, tpat), pick_alt),
      ## .tmat_alt = pmap_chr(list(ref, .alts_happy, tmat), pick_alt),
      .q1_alt = pmap_chr(list(ref, .alts_happy, q1), pick_alt),
      .q2_alt = pmap_chr(list(ref, .alts_happy, q2), pick_alt),
      .dpat_alt = pmap_chr(list(ref, .alts, dpat), pick_alt),
      .dmat_alt = pmap_chr(list(ref, .alts, dmat), pick_alt),
      .trimlen = str_length(ref) - str_length(ref_happy),
      .ref = str_sub(ref, 0, -1 - .trimlen),
      .dpat_alt_trimmed = str_sub(.dpat_alt, 0, -1 - .trimlen),
      .dmat_alt_trimmed = str_sub(.dmat_alt, 0, -1 - .trimlen),
      .alt = pmap_chr(
        list(
          .category, alt_happy, .alts, .alts_happy, .dpat_alt, .dmat_alt,
          .dpat_alt_trimmed, .dmat_alt_trimmed, .trimlen
        ),
        fix_alt
      ),
      # only used for extra_het case
      .extra_is_paternal = case_when(
        .q2_alt == .dpat_alt_trimmed ~ FALSE,
        .q2_alt == .dmat_alt_trimmed ~ TRUE,
        TRUE ~ NA
      ),
      .extra_alt = if_else(
        .category == "extra_het",
        if_else(.extra_is_paternal, .dpat_alt, .dmat_alt),
        NA
      ),
      .query_gt = case_when(
        .category == "extra_het" ~ if_else(.extra_is_paternal, "0|1", "1|0"),
        .category == "extra_null" ~ ".",
        TRUE ~ NA
      ),
      .extra_struct = pmap(
        list(.category, start, ref, .extra_alt, .extra_is_paternal,
             .dpat_alt, .dmat_alt, dpat, dmat),
        build_extra_struct
      ),
      .missing_struct = pmap(
        list(start, ref_happy, .alts_happy, .allow2),
        build_missing_struct
      )
    ) %>%
    # NOTE order matters...ALOT (see fix_var_group)
    arrange(chrom, .vargroup, .category)

  # Loop through each variant group and attempt to fix
  df_fixed <- df_precomputed %>%
    group_by(chrom, .vargroup) %>%
    group_map(fix_var_group) %>%
    bind_rows() %>%
    mutate(.success = .pat_balanced & .mat_balanced & .all_matched & .all_phased)

  # Test that we did our job correctly
  df_fixed %>%
    filter(!.success) %>%
    ensure_empty("Some variant groups could not be solved")
  df_fixed %>%
    filter(!(.ref == ref_happy & .alt == alt_happy)) %>%
    ensure_empty("Some variant groups do not have matching ref and alt")

  df_fixed %>%
    mutate(
      ref = .ref,
      alt = .alt,
      query_gt = .query_gt
    ) %>%
    select(-starts_with("."))
}

# Fix the remaining alleles, which are assumed to be in "groups" that need to
# be handled in a relational manner
fix_groups <- function(df) {
  # Begin by grouping the incoming dataframe kinda like bedtools intersect.
  # All variants with overlapping happy ref alleles will be part of the same
  # group
  df_grouped <- df %>%
    mutate(
      # fix a few remaining star alleles which happen to coincide with null (".")
      # query GTs.
      .nullstar = str_detect(alt, "\\*$") & query_gt == ".",
      dipcall_gt = case_when(
        dipcall_gt == "1|2" & .nullstar ~ "1|0",
        dipcall_gt == "2|1" & .nullstar ~ "0|1",
        TRUE ~ dipcall_gt
      ),
      dpat = as.integer(str_extract(dipcall_gt, "([0-9])(\\/|\\|)([0-9])", 1)),
      dmat = as.integer(str_extract(dipcall_gt, "([0-9])(\\/|\\|)([0-9])", 3)),
      alt = if_else(.nullstar, str_replace(alt, ",\\*$", ""), alt),
      # count the difference between the number of alleles that should be on a
      # given line vs the number of alleles encoded by dipcall's GT.
      .nallele = case_when(
        query_gt == "0/1" & dipcall_gt == "./." ~ -1,
        query_gt == "0/1" ~ 1,
        query_gt == "2/1" ~ -2,
        query_gt == "." & dipcall_gt %in% c("0|1", "0|2", "1|0", "2|0") ~ 1,
        query_gt == "." ~ 2
      )
    ) %>%
    select(-.nullstar) %>%
    arrange(chrom, start, desc(dipcall_gt)) %>%
    mutate(.end = str_length(ref_happy) + start) %>%
    group_by(chrom) %>%
    mutate(.vargroup = cumsum(start > lag(cummax(.end), default = TRUE))) %>%
    group_by(chrom, .vargroup) %>%
    mutate(.groupsize = n()) %>%
    relocate(chrom, start, .end, .vargroup, .groupsize) %>%
    ungroup()

  # Some pairs of mismatched dipcall and query GTs are actually homs which are
  # very simple to fix by simply removing the alt refered to by the dipcall
  # GT assuming it isn't refered to by the truth. This obviously assumes the
  # het being removed is in the other variant of the pair (which we check)
  df_hom_pairs <- df_grouped %>%
    filter(.groupsize == 2) %>%
    filter(query_gt == "1/1" | dipcall_gt %in% c("1|1", "2|2")) %>%
    mutate(
      .alts = str_split(alt, ","),
      .dip_alt = pmap_chr(list(ref, .alts, dpat), pick_alt),
      .alt = pmap_chr(
        list(.alts, pmax(tpat, tmat), dpat),
        function(a, t, d) {
          .a <- ifelse(t < d, head(a, d - 1), a)
          str_flatten(.a, ",")
        }),
      # sanity check, this isn't needed to compute the final alt but if all the
      # assumptions going into this hold, then the shifted dipcall alt should
      # be the same as the existing query alt
      .dip_alt_shift = lag(.dip_alt),
      .query_alt = pmap_chr(list(ref, .alts, q1), pick_alt),
      alt = if_else(is.na(dpat), alt, .alt),
      query_gt = if_else(query_gt == ".", ".", "1|1")
    ) 

  # test the above assumption
  df_hom_pairs %>%
    filter(
      !(
        (.dip_alt_shift == .query_alt | (is.na(.dip_alt_shift) & is.na(.query_alt))) &
          (ref == ref_happy)
      )
    ) %>%
    ensure_empty("unpaired hom variant detected")

  # Anything with groupsize 1 that has a hom happy GT and a null dipcall GT is
  # assumed to have a corresponding dipcall variant that was "missed" in the
  # merge. Phase these cases trivially
  df_mono_hom <- df_grouped %>%
    filter(query_gt == "1/1" & dipcall_gt == "./." & .groupsize == 1) %>%
    mutate(query_gt = "1|1")

  # These are not fixable as the alleles are not balanced
  df_unfixable <- df_grouped %>%
    group_by(chrom, .vargroup) %>%
    mutate(total = sum(.nallele)) %>%
    filter(total != 0) %>%
    ungroup()

  # Fix all groups which are balanced (all of these should succeed)
  df_balanced <- df_grouped %>%
    anti_join(df_hom_pairs, by = "id") %>%
    anti_join(df_mono_hom, by = "id") %>%
    anti_join(df_unfixable, by = "id") %>%
    fix_balanced_groups()

  list(
    fixed = bind_rows(
      df_balanced,
      df_mono_hom,
      df_hom_pairs
    ),
    unfixable = df_unfixable
  )
}

fix_all <- function(df) {
  hom_result <- fix_matching_homs(df)
  het_result <- fix_matching_hets(hom_result$remainder)
  star_result <- fix_star_alleles(het_result$remainder)
  grouped_result <- fix_groups(star_result$remainder)
  list(
    fixed = bind_rows(
      hom_result$final,
      het_result$final,
      star_result$final,
      grouped_result$fixed
    ),
    unfixable = grouped_result$unfixable
  )
}

write_vcf <- function(df, meta, out) {
  df %>%
    arrange(id) %>%
    left_join(meta, by = "id") %>%
    select(chrom, start, id, ref_happy, alt_happy, qual, filter, info, format, truth, query) %>%
    write_tsv(out, col_names = FALSE)
}

read_result <- read_data("small_happy.vcf.gz", "merged_hprc_gt.vcf.gz")
fix_result <- fix_all(read_result$gt)

ensure(
  nrow(fix_result$fixed) + nrow(fix_result$unfixable) == nrow(read_result$gt),
  "final and original row number mismatch"
)
    
fix_result$fixed %>%
  mutate(
    query = str_replace(query, "^[^:]+", ""),
    query = paste0(query_gt, query)
  ) %>%
  write_vcf(read_result$meta, "fixed.vcf.gz")

fix_result$unfixable %>%
  write_vcf(read_result$meta, "notfixed.vcf.gz")
