#' @name mhcI
#' @title Prediction of Nonself Peptide Presentation on Human MHC Class I
#' @description It determines which non-self peptides can be presented by a given HLA class I allele. This function takes a sequence for a stimulating antigen and the corresponding self antigen, and given a defined sequence length, queries the IEDB API with the user's choice of peptide binding prediction method. The set of peptides present in the results for the stimulating antigen but not the self antigen are then carried forward as non-self peptides. If desired, the user can adjust the default thresholds (by IC50 binding affinity or percentile rank) used to define "strong" and "weak" binders. The output is a dataframe of non-self peptides that are predicted to bind to the presenting allele.
#' @param ag_present
#' character vector, presenting allele, formatting examples - A*01:01 or HLA-A*01:01
#' @param ag_stim
#' character vector, stimulating antigen, can either be an HLA class I allele entered in the same format as ag_present, or a character vector of the amino acid sequence of the protein
#' @param ag_self
#' character, self antigen, can either be an HLA class I allele entered in the same format as ag_present, or a character vector of the amino acid sequence of the protein
#' @param seq_len
#' string, length of peptides to consider
#' @param fd_out
#' string, output folder name, set to tempdir()
#' @param method
#' string, IEDB prediction method to be used. Options are ann, comblib_sidney2008, consensus, netmhccons, netmhcpan, netmhcpan_ba
#' netmhcpan_el, netmhcstabpan, pickpocket, recommended, smm, smmpmbec. Default is netmhcstabpan
#' @param thold_ic50
#' vector. Defines the thresholds required to be included in results, and to be labeled, "strong" or "weak" binder. Multiple prediction methods are used, each of which provide different raw outputs (i.e. IC50, "strength", "score"). Our justification for the default thresholds is listed in the mhcI_hu vignette, however the user may choose to specify alternate cutoffs if desired.
#' @param thold_pect_rank
#' vector, IEDB adjusts all outputs in comparison to a set of random natural peptides in order to determine an normalized adjusted percentile rank. With normalized ranks, the same thresholds can be used across different methods. Default thresholds are <2\% for strong binders and <10\% for weak binders.
#' @param url_iedb
#' string, iedb api url
#' @param noneself
#' logic, whether to subtract self peptide from stim. If TRUE|T, self peptides will be subtracted from stim. If FALSE|F no self subtraction.
#'
#' @export
#' @return data frame, MHC I binding prediction result table
#' @import
#' dplyr
#' fs
#' httr
#' janitor
#' stringr
#' tibble
#' tidyverse
#' @importFrom
#' stats setNames
#' @importFrom
#' rlang as_string
#' @examples
#' \donttest{
#' mhcI(ag_present=c("HLA-A*03:01"),
#' ag_stim=c("A_01_01","A_02_06"),
#' ag_self=c("B_07_03"),
#' method = "rec")
#' }

mhcI <- function(ag_present,
                 ag_stim,
                 ag_self = "",
                 seq_len = '9',
                 fd_out = as.character(paste0(tempdir(), "/", "outputs", "/")),
                 method = "consensus",
                 thold_ic50 = c(50, 500),
                 thold_pect_rank = c(2, 10),
                 url_iedb = 'http://tools-cluster-interface.iedb.org/tools_api/mhci/',
                 noneself = TRUE) {
  #* pre check: method and folder *#
  # allow ambiguous string of method
  if (str_detect(method, "rec")) {
    method <- "recommended"
  } else if (str_detect(method, "cons|con")) {
    method <- "consensus"
  } else if (str_detect(method, "comblib")) {
    method <- "comblib_sidney2008"
  } else if (str_detect(method, "_ba|ba")){
    method <- "netmhcpan_ba"
  } else if (str_detect(method, "_el|el")){
    method <- "netmhcpan_el"
  } else {
    method <- "consensus"
  }

  if (!requireNamespace("httr", quietly = TRUE)) {
    stop(
      "Package \"httr\" must be installed to use this function.",
      call. = FALSE
    )
  }
  #* end of pre check *#

  #* 1. setup: ref table and folders *#
  ref_hu <- read.table(gzfile(system.file("extdata/ref", "ref_human_mhcI.gz", package = "epitopeR")), row.names=1)

  # make sure noneself is false if ag_self is null
  if (ag_self == "" | is.na(ag_self)) {
    noneself <- FALSE
  }

  # if output folder doesn't exist, create one; else, clean up the folder
  if(!dir.exists(fd_out)) {
    dir.create(file.path(fd_out))
  } else {
    file.remove(fs::dir_ls(fd_out, glob = "*.txt"))
  }
  #* end of setup *#

  #* 2. validation and preprocess  *#
  # validate name of presenting antigen
  val_ag_name(ag_present)

  # this evaluation only apply to mhcI
  if (length(ag_present) != length(seq_len)) {
    warning("The same number of alleles, lengths, must be provided.\n")
    warning("Prediction is going to made only on ", ag_present[1] , " with length ", seq_len[1], "\n")
    ag_present <- ag_present[1]
    seq_len <- seq_len[1]
  }

  # if input stim/self are sequence text
  if(all(unique(c(unlist(strsplit(toupper(ag_self), split = "")),
                  unlist(strsplit(toupper(ag_stim), split = "")))) %in%
         c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))) {

    if(ag_stim == ag_self) {
      warning("stim antigens are identical with self antigens, no prediction will be made.\n")
    }

  } else {

    # if human allele names
    ag_self <- preproc(allele_in = ag_self, link = url_iedb)
    # keep the elements of stim that are not in self for computational efficiency
    # this is different concept with find_nonself()
    ag_stim <- preproc(allele_in = ag_stim, link = url_iedb)
    ag_stim <- setdiff(ag_stim, ag_self)

    if(all(ag_stim %in% ag_self)) {
      warning("stim antigens are identical with self antigens, no prediction will be made.\n")
    }

    if(all(ag_present %in% ag_stim)) {
      warning("presenting antigens are identical with stim antigens, no prediction will be made.\n")
    }
  }
  #* end of preprocess *#

  #* 3. prep for api calls *#
  if (!all(str_detect(c(ag_self, ag_stim), "[0-9]"))) {
    tmp_self <- pull_obj_name(ag_self) %>% as_string()
    tmp_stim <- pull_obj_name(ag_stim) %>% as_string()
    ref_self_stim <- as.data.frame(c(ag_self, ag_stim)) %>%
      mutate(allele = c(tmp_self,
                        tmp_stim),
            iedb_query_seq = c(ag_self, ag_stim),
            allele_query_nm = c(ag_self, ag_stim),
            ag_type = c(rep("self", length(ag_self)),
                        rep("stim", length(ag_stim)))) %>%
      select(allele, iedb_query_seq, allele_query_nm, ag_type)
  } else {
    ref_self_stim <- pull_seq(alleles_in = c(ag_self, ag_stim),
                              tbl_ref_in = ref_hu) %>%
      arrange(factor(allele, levels = c(ag_self, ag_stim))) %>% # !important: need reorder rows by input antigens
      mutate(iedb_query_seq = seq,
             allele_query_nm = gsub("_", ":", sub("_", "*", allele)),
             ag_type = c(rep("self", length(ag_self)),
                         rep("stim", length(ag_stim)))) %>%
      select(allele, iedb_query_seq, allele_query_nm, ag_type)
  }

  seq_names <- ref_self_stim %>%
    rowid_to_column() %>%
    mutate(antigen = allele,
           length = nchar(iedb_query_seq),
           allele = str_replace_all(antigen, "_", "*"),
           seq_num = rowid) %>%
    select(seq_num, allele, length, antigen, iedb_query_seq, ag_type)

  seq_names_sht <- seq_names %>%
    select(antigen, ag_type)

  #* end of prep for api calls *#

  #* 4. iedb prediction *#
  for(j in 1:nrow(seq_names)) { # for each antigen
    # set parameter list
    param2httr = list('method' = method,
                      'sequence_text' = as.character(seq_names$iedb_query_seq[j]),
                      'allele' = str_c(ag_present, collapse = ","),
                      'length' = seq_len)
    # send request to api
    res_api <- httr::POST(url = url_iedb,
                          body = param2httr)

    # extract content from the request
    tmp_nm <- ifelse(nchar(seq_names$antigen[j]) > 20, str_sub(seq_names$antigen[j], 1,20),
                     seq_names$antigen[j])
    cat(httr::content(res_api, "text"),
        file = paste0(fd_out, tmp_nm, "_", method, ".txt"))
    # cat(httr::content(res_api, "text"),
    #     file = paste0(fd_out, seq_names$antigen[j], "_", method, ".txt"))
  }
  rm(j)
  #* end of prediction *#

  #* 5. combine all prediction result into one *#
  tmp <- comb_pred_tbl_mhcI(nm_method = method,
                            nm_fd = fd_out,
                            thold_score = thold_ic50,
                            thold_rank = thold_pect_rank)
  #* end of combine *#

  #* 6. peptide table *#
  if (dim(tmp)[1] == 0) {
    final <- data.frame()
  } else {
    tmp <- left_join(tmp, seq_names_sht, by = "antigen")
    if (noneself) {
      final <- tmp %>%
        find_nonself()
      if (all(length(ag_stim) > 0, length(ag_self) > 0,  dim(final)[1] == 0)) {
        message("No non-self peptides detected.")
      }
    } else {
      final <- tmp
    }
  }
  #* end of peptide table *#

  return(final)
}
