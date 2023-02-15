#' Basic functions of the package
#'
#' val_ag_name() validation of name of binding locus
#' @param ag_present a vector of locus name(s) to make binding predictions for
#'
#' preproc() format and validate allele names
#' @param allele_in a vector contains allele name(s)
#' @param link a string of url of MHC I or II api
#'
#' pull_seq() pull out sequence of each allele based on ref table
#' @param alleles_in vector, allele names
#' @param tbl_ref_in dataframe, reference table, default is human_all.csv from github
#'
#' comb_pred_tbl() combine individual prediction tables by method, exclude none-binders and keep strong and weak binders only
#' @param nm_method string, prediction method used for IEDB prediction
#' @param nm_sht string, short name of alleles
#' @param nm_fd string, folder name which contains predict tables from IEDB
#' @param thold_score list of vectors, binder thresholds by score
#' @param thold_rank vector, binder thresholds by rank
#'
#' comb_pred_tbl_mhcI() combine individual human mhcI prediction tables by method, exclude none-binders and keep strong and weak binders only
#' @param nm_method string, prediction method used for IEDB prediction
#' @param nm_fd string, folder name which contains predict tables from IEDB
#' @param thold_score list of vectors, binder thresholds by ic50 score
#' @param thold_rank vector, binder thresholds by percentile rank
#'
#' find_nonself() find nonself binding peptides
#' @param dat_in dataframe, combined prediction tables
#'
#' pull_ag_self() pull out aligned ag_self based on aligned ag_stim position, add a core mutation flag
#' @param vec_in dataframe with pep_stim, core, aligned ag_stim and ag_self columns
#'
#' find_core_mut() find mutation position to core
#' @param dat_in dataframe with pep_stim, core, pep_self selected from pull_ag_self
#'
#' align_seq() align protein sequences
#' @param seq1 string, unaligned sequence of ag_stim
#' @param seq2 string, unaligned sequence of ag_self
#' @param gapopening numeric, the cost for opening a gap in the alignment.
#' @param gapextension numeric, the incremental cost incurred along the length of the gap in the alignment
#'
#' pull_obj_name()
#' @param x, name of an object
#'
#' @name utils
#' @import
#' dplyr
#' fs
#' janitor
#' stringr
#' tidyverse
#' utils
#' readr
#' @importFrom
#' stats setNames
#' @importFrom
#' purrr map_df
#' @importFrom
#' seqinr read.fasta
#' @importFrom
#' Biostrings BString
#' @importFrom
#' Biostrings pairwiseAlignment
#' @importFrom
#' Biostrings AAString
#'
NULL
#> NULL

#' @rdname utils
val_ag_name <- function(ag_present) {
  ag_present <- tolower(ag_present)

  #* start of human mhcII dp dq *#
  ck_dp <- ag_present[str_detect(ag_present, "dpa|apb")]
  ck_dq <- ag_present[str_detect(ag_present, "dqa|dqb")]

  if(length(ck_dp) > 0 & !(all(str_detect(ck_dp, "/")))) {
    warning("Please make sure to input both DPA and DPB alleles separated by /.",
            "\n",
            "ex: DPA1*01:03/DPB1*01:01",
            "\n")
  }

  if(length(ck_dq) > 0 & !(all(str_detect(ck_dq, "/")))) {
    warning("Please make sure to input both DQA and DQB alleles separated by /.",
            "\n",
            "ex: DQA1*01:01/DQB1*02:02",
            "\n")
  }
  #* end of human mhcII dp dq *#

  #* start of human mhcI *#
  ck_mhc1 <- ag_present[str_detect(ag_present, "a\\*|b\\*|c\\*|e\\*|g\\*")]

  if(length(ck_mhc1) > 0 & !(all(str_detect(ck_mhc1, "hla-")))) {
    warning("Please make sure you have input the presenting allele with the prefix HLA-",
            "\n",
            "ex: HLA-A*02:06",
            "\n")
  }
  #* end of human mhcI *#

  #* start of mouse, mhcI and II *#
  ms <- c("H2-IAb", "H2-IAd", "H2-IAk","H2-IAs", "H2-IEd", "H2-IEk", # mouse, mhcII
          "H-2-Db", "H-2-Dd", "H-2-Kb", "H-2-Kd", "H-2-Kk", "H-2-Ld") # mouse, mhcI
  ck_ms <- ag_present[str_detect(ag_present, "h2|h-2")]

  if(length(ck_ms) > 0 & !(all(ck_ms %in% tolower(ms)))) {
    warning("Please choose the presenting allele from ",
            str_c(ms, collapse = ","),
            "\n")
  }
  #* end of mouse *#
}

#' @rdname utils
preproc <- function(allele_in, link) {
  # 1. format conversion: a*01:01 to A_01_01
  allele_out <- toupper(allele_in) %>%
    str_replace_all(., "\\*|\\:", "_") %>%
    unique()

  # 2. letter code of protein sequences
  prt_seq <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

  ck_seq <- unique(unlist(strsplit(allele_out, split = "")))
  ck_seq <- ck_seq[!(ck_seq %in% prt_seq)]

  # 3. validation
  # if input is sequence string but contains invalid character, then STOP
  if(!all(str_detect(allele_out, "[0-9]")) & length(ck_seq) > 0 ) {
    stop(cat("Input sequence contains invalid character. '", ck_seq[1], "'at position", str_locate(allele_out, ck_seq[1])[,1]))
  }

  return(allele_out)
}

#' @rdname utils
pull_seq <- function(alleles_in,
                     tbl_ref_in) {
  tbl_seq <- tbl_ref_in %>% filter(allele %in% c(alleles_in))

  return(tbl_seq)
}

#' @rdname utils
comb_pred_tbl <- function(nm_method, nm_sht, nm_fd, thold_score, thold_rank) {
  # list all of output files, and filter out empty ones (403 Forbidden, file starts with "<!DOCTYPE")
  fl <- dir_ls(nm_fd, glob =  paste0("*",nm_method, ".txt"))
  fl <- fl[sapply(fl, file.size) > 500] %>% as.vector()

  # map all non-empty predict tables into 1
  if (length(fl) == 0) {
    datout <- data.frame()

  } else {
    datout <- fl %>%
      setNames(nm = .) %>%
      map_df(~read_tsv(.x, col_types = cols(), col_names = TRUE), .id = "id") %>%
      as.data.frame()

    tmp_antigen <- data.frame(id = datout$id, str = rep(paste0("_", nm_method, ".txt"), dim(datout)[1]))
    tmp_antigen$antigen <- str_remove(gsub(".*/", "", tmp_antigen$id), tmp_antigen$str)
    datout$antigen <- tmp_antigen$antigen

  ###*** start of netmhciipan ***###
  if (nm_method == "netmhciipan") {
    # score_val: IC50
    datout <- datout %>%
      #mutate(antigen = str_remove(gsub(".*/", "", id), paste0("_", nm_method, ".txt"))) %>%
      select(-c(id, seq_num)) %>%
      filter(ic50 != "-") %>%
      mutate(pep_stim = peptide,
             method = nm_method,
             core = core_peptide,
             score_val = as.numeric(ic50),
             rank_val = as.numeric(rank)) %>%
      mutate(strength_ic50 = case_when(score_val <= min(thold_score$cutoff_netpan) ~ "strong",
                                       score_val <= max(thold_score$cutoff_netpan) ~ "weak",
                                       TRUE ~ "no"),
             strength_rank = case_when(rank_val <= min(thold_rank) ~ "strong",
                                       rank_val <= max(thold_rank) ~ "weak",
                                       TRUE ~ "no")) %>%
      filter(strength_ic50 %in% c("strong", "weak") | strength_rank %in% c("strong", "weak"))

    if (dim(datout)[1] > 0) {
      datout <- datout %>%
        left_join(., nm_sht, by = "antigen") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method) %>%
        distinct()
    } else {
      datout <- data.frame()
    }
  }
  ###*** end of netmhciipan ***###

  ###*** start of recommended ***###
  if (str_detect(nm_method, "rec")) {
    datout <- datout %>%
      #mutate(antigen = str_remove(gsub(".*/", "", id), paste0("_", nm_method, ".txt"))) %>%
      select(-c(id, seq_num)) %>%
      mutate(pep_stim = peptide)

    # start of comblib #
    # score_val: IC50 (column "comblib_score")
    comblib <- datout %>%
      filter(comblib_score != "-")

    if(dim(comblib)[1] > 0 ) {
      comblib <- comblib %>%
        mutate(method = "comblib",
               core = comblib_core,
               score_val = as.numeric(comblib_score),
               rank_val = as.numeric(comblib_rank)) %>%
        mutate(strength_ic50 = case_when(score_val <= min(thold_score$cutoff_comblib) ~ "strong",
                                         score_val <= max(thold_score$cutoff_comblib) ~ "weak",
                                         TRUE ~ "no"),
               strength_rank = case_when(rank_val <= min(thold_rank) ~ "strong",
                                         rank_val <= max(thold_rank) ~ "weak",
                                         TRUE ~ "no")) %>%
        filter(strength_ic50 %in% c("strong", "weak") | strength_rank %in% c("strong", "weak")) %>%
        left_join(., nm_sht, by = "antigen") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method)  %>%
        distinct()
    } else {
       comblib <- data.frame()
      }
    # end of comblib #

    # start of nn_align #
    # score_val: IC50
    nn_align <- datout %>%
      filter(nn_align_ic50 != "-")

    if (dim(nn_align)[1] > 0) {
      nn_align <- nn_align %>%
        mutate( method = "nn_align",
                core = nn_align_core,
                score_val = as.numeric(nn_align_ic50),
                rank_val = as.numeric(nn_align_rank)) %>%
        mutate(strength_ic50 = case_when(score_val <= min(thold_score$cutoff_nn_align) ~ "strong",
                                         score_val <= max(thold_score$cutoff_nn_align) ~ "weak",
                                         TRUE ~ "no"),
               strength_rank = case_when(rank_val <= min(thold_rank) ~ "strong",
                                         rank_val <= max(thold_rank) ~ "weak",
                                         TRUE ~ "no")) %>%
        filter(strength_ic50 %in% c("strong", "weak") | strength_rank %in% c("strong", "weak")) %>%
        left_join(., nm_sht, by = "antigen") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method) %>%
        distinct()
    } else {
       nn_align <- data.frame()
      }
    # end of nn_align #

    # start of smm_align #
    # score_val: IC50
    smm_align <- datout %>%
      filter(smm_align_ic50 != "-")

    if (dim(smm_align)[1] > 0) {
      smm_align <- smm_align %>%
        mutate( method = "smm_align",
                core = smm_align_core,
                score_val = as.numeric(smm_align_ic50),
                rank_val = as.numeric(smm_align_rank)) %>%
        mutate(strength_ic50 = case_when(score_val <= min(thold_score$cutoff_smm_align) ~ "strong",
                                         score_val <= max(thold_score$cutoff_smm_align) ~ "weak",
                                         TRUE ~ "no"),
               strength_rank = case_when(rank_val <= min(thold_rank) ~ "strong",
                                         rank_val <= max(thold_rank) ~ "weak",
                                         TRUE ~ "no")) %>%
        filter(strength_ic50 %in% c("strong", "weak") | strength_rank %in% c("strong", "weak")) %>%
        left_join(., nm_sht, by = "antigen") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method) %>%
        distinct()
    } else {
       smm_align <- data.frame()
    }
    # end of smm_align #

    # start of sturniolo #
    # score_val: sturniolo_score; only label strong binders
    sturniolo <- datout %>%
      filter(sturniolo_score != "-")

    if (dim(sturniolo)[1] > 0) {
      sturniolo <- sturniolo %>%
        mutate(method = "sturniolo",
               core = sturniolo_core,
               score_val = as.numeric(sturniolo_score),
               rank_val = as.numeric(sturniolo_rank)) %>%
        mutate(strength_ic50 = case_when(score_val >= max(thold_score$cutoff_sturniolo) ~ "strong",
                                         TRUE ~ "no"),
               strength_rank = case_when(rank_val <= min(thold_rank) ~ "strong",
                                         rank_val <= max(thold_rank) ~ "weak",
                                         TRUE ~ "no"))  %>%
        filter(strength_ic50 %in% c("strong") | strength_rank %in% c("strong")) %>%
        left_join(., nm_sht, by = "antigen") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method) %>%
        distinct()
    } else {
       sturniolo <- data.frame()
      }
    # end of sturniolo #

     if(dim(comblib)[1] > 0 | dim(nn_align)[1] > 0 | dim(smm_align)[1] > 0 | dim(sturniolo)[1] > 0) {
      datout <- rbind(comblib, nn_align, smm_align, sturniolo) %>%
         group_by(pep_stim) %>%
         filter(rank_val == min(rank_val)) %>%
         ungroup()
     } else {
       datout <- data.frame()
    }

  }
  ###*** end of recommended ***###

  ###*** start of netpna el ***###
  if (str_detect(nm_method, "el")) {
    # score_val: percentile rank
    datout <- datout %>%
      #mutate(antigen = str_remove(gsub(".*/", "", id), paste0("_", nm_method, ".txt"))) %>%
      select(-c(id, seq_num)) %>%
      filter(rank != "-") %>%
      mutate(pep_stim = peptide,
             method = nm_method,
             core = core_peptide,
             score_val = as.numeric(score),
             rank_val = as.numeric(rank)) %>%
      mutate(strength_ic50 = "na",
             strength_rank = case_when(rank_val <= min(thold_rank) ~ "strong",
                                       rank_val <= max(thold_rank) ~ "weak",
                                       TRUE ~ "no")) %>%
      filter(strength_rank %in% c("strong", "weak"))

    if (dim(datout)[1] > 0) {
      datout <- datout %>%
        left_join(., nm_sht, by = "antigen") %>%
        select(allele, antigen, ag_type, start, end, pep_stim, core, length, strength_ic50, score_val, strength_rank, rank_val, method) %>%

        distinct()
    } else {
      datout <- data.frame()
    }
  }
  ###*** end of netpan el ***###

  }

  return(datout)
}

#' @rdname utils
comb_pred_tbl_mhcI <- function(nm_method, nm_fd, thold_score, thold_rank) {
  #* step 1. map all prediction tables *#
  # if no valid prediction table, then return empty data frame; else continue
  fl <- dir_ls(nm_fd, glob =  paste0("*", nm_method, ".txt"))
  fl <- fl[sapply(fl, file.size) > 500] %>% as.vector()

  # if no file, the data out is empty, else -
  #   if consensus, then output from wide to long; else -
  #      work on ic50 or score/percentile_rank

  #* start of outer if *#
  if (length(fl) == 0 ) {
    datout <- data.frame()
    warning("No predictions were made. Please validate allele names and length.")
  } #* end of outer if *#
  else { #* start of big outer else *#
    datout <- fl %>%
      setNames(nm = .) %>%
      map_df(~read_tsv(.x, col_types = cols(), col_names = TRUE), .id = "id") %>%
      rename_with(str_to_lower) %>%
      mutate(antigen = str_remove(str_remove(gsub(".*/", "", id), paste0("_", nm_method, ".txt")), "_netmhcpan")) %>%
      select(-id, -seq_num) %>%
      as.data.frame()

    #* start of if consensus *#
    if (nm_method == "consensus") {
      #* start of ann *#
      ann <- datout %>%
        filter(ann_ic50 != "-")

      if(dim(ann)[1] > 0) {
        ann <- ann %>%
          mutate(method = "ann",
                 ic50 = ann_ic50,
                 percentile_rank = "",
                 pep_stim = peptide) %>% # set a fake score value
          mutate(strength_ic50 = case_when(ic50 <= min(thold_score) ~ "strong",
                                           ic50 <= max(thold_score) ~ "weak",
                                           TRUE ~ "no"),
                 strength_rank = "") %>%
          filter(strength_ic50 %in% c("strong", "weak")) %>%
          select(allele, start, end, length, pep_stim, antigen,
                 ic50, percentile_rank, strength_ic50, strength_rank, method) %>%
          distinct()
      } else {
        ann <- data.frame()
      }
      #* end of ann *#

      #* start of smm *#
      smm <- datout %>%
        filter(smm_ic50 != "-")

      if(dim(smm)[1] >0) {
        smm <- smm %>%
          mutate(method = "smm",
                 ic50 = smm_ic50,
                 percentile_rank = "",
                 pep_stim = peptide) %>% # set a fake score value
          mutate(strength_ic50 = case_when(ic50 <= min(thold_score) ~ "strong",
                                           ic50 <= max(thold_score) ~ "weak",
                                           TRUE ~ "no"),
                 strength_rank = "") %>%
          filter(strength_ic50 %in% c("strong", "weak")) %>%
          select(allele, start, end, length, pep_stim,  antigen,
                 ic50, percentile_rank, strength_ic50, strength_rank, method) %>%
          distinct()
      } else{
        smm <- data.frame()
      }
      #* end of smm *#

      #* start of comblib *#
      comblib <- datout %>%
        filter(comblib_sidney2008_score != "-")

      if(dim(comblib)[1] > 0){
        comblib <- comblib %>%
          mutate(method = "comblib",
                 percentile_rank = comblib_sidney2008_rank,
                 ic50 = "",
                 pep_stim = peptide) %>% # set a fake ic50 value
          mutate(strength_ic50 = "",
                 strength_rank = case_when(percentile_rank <= min(thold_rank) ~ "strong",
                                           percentile_rank <= max(thold_rank) ~ "weak",
                                           TRUE ~ "no")) %>%
          filter(strength_rank %in% c("strong", "weak")) %>%
          select(allele, start, end, length, pep_stim,  antigen,
                 ic50, percentile_rank, strength_ic50, strength_rank, method) %>%
          distinct()
      } else {
        comblib <- data.frame()
      }

      if((dim(ann)[1] > 0 | dim(smm)[1] > 0) & dim(comblib)[1] > 0) {
        tmp1 <- rbind(ann, smm) %>%
          group_by(pep_stim) %>%
          filter(ic50 == min(ic50)) %>%
          ungroup()

        tmp2 <- comblib %>% filter(!(pep_stim %in% tmp1$pep_stim))

        datout <- rbind(tmp1, tmp2)
      } else if((dim(ann)[1] > 0 | dim(smm)[1] > 0) & dim(comblib)[1] == 0) {
        datout <- rbind(ann, smm) %>%
          group_by(pep_stim) %>%
          filter(ic50 == min(ic50)) %>%
          ungroup()
      } else if((dim(ann)[1] == 0 & dim(smm)[1] == 0) & dim(comblib)[1] > 0) {
        datout <- comblib
      }
      else {
        datout <- data.frame()
      }
    } #* end of if consensus *#

    #* start of inner else *#
    else {
      # if output contains ic50
      if (dim(datout)[1] > 0 & 'ic50' %in% colnames(datout)) {
        datout <- datout %>%
          filter(ic50 != "-") %>%
          mutate(pep_stim = peptide,
                 method = nm_method) %>%
          mutate(strength_ic50 = case_when(ic50 <= min(thold_score) ~ "strong",
                                           ic50 <= max(thold_score) ~ "weak",
                                           TRUE ~ "no")) %>%
          filter(strength_ic50 %in% c("strong", "weak")) %>%
          select(allele, start, end, length, pep_stim,  antigen,
                 ic50, percentile_rank, method)
      }

      # use percentile rank if output doesn't contain ic50
      if (dim(datout)[1] > 0 & !('ic50' %in% colnames(datout)) & 'score' %in% colnames(datout)) {
        datout <- datout %>%
          filter(score != "-") %>%
          mutate(pep_stim = peptide,
                 method = nm_method) %>%
          mutate(strength_rank = case_when(percentile_rank <= min(thold_rank) ~ "strong",
                                           percentile_rank <= max(thold_rank) ~ "weak",
                                           TRUE ~ "no")) %>%
          filter(strength_rank %in% c("strong", "weak"))
      }
    }  #* start of inner else *#

  }#* end of outer else *#

  return(datout)
}

#' @rdname utils
find_nonself <- function(dat_in) {
    self_peps <- dat_in %>%
                  filter(ag_type == "self")
    allo_peps <- dat_in %>%
                  filter(ag_type == "stim") %>%
                  anti_join(self_peps, by = "pep_stim") %>%
      select(-ag_type)
    return(allo_peps)
}

#' @rdname utils
pull_ag_self <- function(vec_in) {
  nmer <- vec_in$length
  # for each pep_stim, find starting position of best matches to aligned ag_stim, max allowed mismatch is 10
  vec_out <- Biostrings::matchPattern(vec_in$pep_stim,
                                      BString(vec_in$ag_stim),
                                      fixed=TRUE,
                                      max.mismatch = 10)@ranges %>%
    data.frame()

  if(dim(vec_out)[1] > 0) {
    vec_out <- vec_out %>%
      rowwise() %>%
      # for each range of best match result
      mutate(pep_stim = vec_in$pep_stim,
             ag_stim = vec_in$ag_stim,
             ag_self = vec_in$ag_self,
             core = vec_in$core) %>%
      mutate(pep1 = str_sub(ag_stim, start, end), # substr based on starting position
             cnt = str_count(pep1, "-")) %>% # then count number of dashes of each substr
      # remove existing dashes and fill out removed ones with number of none-dashes chars
      # we only do this once to simplify the process
      mutate(pep2 = str_sub(ag_stim, start, end+cnt), "-") %>%
      # find distance between each pair of substr and pep_stim, and keep the closest one
      mutate(cnt_match = adist(pep_stim, pep2)) %>%
      data.frame() %>%
      filter(cnt_match == min(cnt_match)) %>%
      # pull out substr of pep_stim and pep_self, it could be w or w/o dashes
      mutate(tmp_cnt = start + as.numeric(nmer) - 1) %>%
      mutate(pep_stim_alg = str_sub(ag_stim, start, tmp_cnt),
             pep_self = str_sub(ag_self, start, tmp_cnt)) %>%
      # flag of yes or no of dash
      mutate(dash = ifelse(str_detect(pep_stim_alg, "-") | str_detect(pep_self, "-"),
                           "yes",
                           "no")) %>%
      mutate(start = vec_in$start,
             end = vec_in$end) %>%
      select(pep_stim, pep_self, core, start, end, dash)
  } else {
    vec_out <- vec_in %>%
      mutate(pep_stim = vec_in$pep_stim,
             pep_self = "",
             dash = "not found") %>%
      select(pep_stim, pep_self, core, start, end, dash)
  }

  return(vec_out)
}

#' @rdname utils
find_core_mut <- function(dat_in) {
  # step 1. constants
  # length and number of peptide
  nmer <- unique(nchar(dat_in$pep_stim))
  num_pep <- length(dat_in$pep_stim)

  # step 2: add start and end positions of core to pep_stim
  dat_in <- dat_in %>%
    rowwise() %>%
    mutate(core_st = str_locate(pep_stim, core)[1],
           core_ed = str_locate(pep_stim, core)[2])

  dat_na <- dat_in %>%
    filter(is.na(core_st)) %>%
    mutate(core_mut = "not found") %>%
    select(pep_stim, pep_self, core, start, end, core_mut)

  dat_in <- dat_in %>% filter(!is.na(core_st))

  if(dim(dat_in)[1] > 0) {
    # step 3: split pep_stim and pep_self into char, and combinde them into one matrix
    tmp_stim <- str_split_fixed(dat_in$pep_stim, "", nmer) %>%
      data.frame() %>%
      setNames(c(paste0("stim", seq(1:nmer))))

    tmp_self <- str_split_fixed(dat_in$pep_self, "", nmer) %>%
      data.frame() %>%
      setNames(c(paste0("self", seq(1:nmer))))

    mutat <- cbind(tmp_stim, tmp_self)

    # step 4: for each pep_stim, find mutation position
    for(i in 1:num_pep){
      #for (j in 1:length(nmer)){ # this line is a bug :-)
      for (j in 1:nmer){
        if(mutat[i,j] != mutat[i, j+nmer]){
          #mutat[i,j] = j # it's a consequential bug
          mutat[i,j] = 1
        } else {
          mutat[i,j] = 0
        }
      }
    }

    mutat <- mutat %>%
      select(names(.)[str_detect(names(.), "stim")]) %>%
      mutate_all(as.numeric)

    # step 5: add flag of mutation to core position
    dat_out <- cbind(dat_in, mutat) %>% mutate(core_mut = 0)

    for(i in 1:num_pep) {
      tmp_nm <- paste0("stim",seq(dat_out[i, ]$core_st[1], dat_out[i, ]$core_ed[1], 1))
      tmp <- dat_out[i, ] %>% select(any_of(tmp_nm)) %>% mutate(core_mut = rowSums(.)) %>% select(core_mut)

      dat_out[i,]$core_mut <- tmp$core_mut
    }

    dat_out <- dat_out %>%
      rowwise() %>%
      mutate(core_mut = ifelse(core_mut > 0, "yes", "no")) %>%
      select(pep_stim, pep_self, core, start, end, core_mut)
  }
  if(dim(dat_in)[1] > 0 & dim(dat_na)[1] > 0){
    final <- rbind(dat_out, dat_na)
  }

  if(dim(dat_in)[1] > 0 & dim(dat_na)[1] == 0){
    final <- dat_out
  }

  if(dim(dat_in)[1] == 0 & dim(dat_na)[1] > 0){
    final <- dat_na
  }

  return(final)
}

##' @rdname utils
align_seq <- function(seq1, seq2, gapopening = 0, gapextension = 8) {

  algn <- pairwiseAlignment(pattern = AAString(seq2),
                            subject = AAString(seq1),
                            type = "local", # align string fragments
                            substitutionMatrix = "BLOSUM50",
                            gapOpening = gapopening,
                            gapExtension = gapextension)

  seq1 <- toupper(toString(algn@subject))
  seq2 <- toupper(toString(algn@pattern))

  return(c(seq1_algn = seq1, seq2_algn = seq2))
}

##' @rdname utils
pull_obj_name <- function(x) {
  #* the function is modified based on Inner() on stackoverflow, authored by BrodieG *#
  #* step 1: label and evaluate the expression *#
  # label the expression, quote the expression for re-use in the later step
  my.call <- quote(substitute(x))

  # evaluate the expression, and get the actual object
  var.name <- eval(my.call)

  #* step 2: reversibly sys.frames calls *#
  for(i in rev(head(sys.frames(), -1L))) {
    my.call[[2]] <- var.name         # this is where we re-use it, modified to replace the variable
    var.name <- eval(my.call, i)
  }

  #* last step: cast symbol to character string, then return *#
  #return(as_string(var.name))
  return(var.name)
}

