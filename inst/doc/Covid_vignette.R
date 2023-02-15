## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  eval = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----lib, eval = TRUE---------------------------------------------------------
suppressPackageStartupMessages({
library(epitopeR)
library(tidyverse)
library(ggVennDiagram)
})

## ----spike, eval = TRUE-------------------------------------------------------
spike <- "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"

## ----mhc1, eval = TRUE--------------------------------------------------------
spike_peps <- mhcI(ag_present = "HLA-A*02:01", ag_stim = spike)

## ----a2 pep, eval = TRUE------------------------------------------------------
a2_spike_peps <- c("YLQPRTFLL", "RLQSLQTYV","VLNDILSRL",  "RLNEVAKNL", "NLNESLIDL", "FIAGLIAIV", "TLDSKTQSL", "VVFLHVTYV", "KLPDDFTGCV", "SIIAYTMSL", "KIADYNYKL", "LLFNKVTLA","ALNTLVKQL", "RLDKVEAEV",  "LITGRLQSL",  "RLITGRLQSL", "FLHVTYVPA",  "KIYSKHTPI")

## ----comp, eval = TRUE, out.width='100%', dpi = 300---------------------------
spike_peps_vector <- as.vector(unique(spike_peps$pep_stim))

ggVennDiagram(x = list(spike_peps_vector, a2_spike_peps), 
              category.names = c("Predicted Peptides", "Known Peptides"), 
              set_size = 2,
              label = "count") + 
              theme(legend.position = "none")


## ----omicron spike, eval = TRUE-----------------------------------------------
omicron_spike <- spike

str_sub(omicron_spike, 67, 67) <- "V"
str_sub(omicron_spike, 95, 95) <- "I"
str_sub(omicron_spike, 145, 145) <- "D"
str_sub(omicron_spike, 212, 212) <- "I"
str_sub(omicron_spike, 339, 339) <- "D"
str_sub(omicron_spike, 371, 371) <- "L"
str_sub(omicron_spike, 373, 373) <- "P"
str_sub(omicron_spike, 375, 375) <- "F"
str_sub(omicron_spike, 417, 417) <- "N"
str_sub(omicron_spike, 440, 440) <- "K"
str_sub(omicron_spike, 446, 446) <- "S"
str_sub(omicron_spike, 477, 477) <- "N"
str_sub(omicron_spike, 478, 478) <- "K"
str_sub(omicron_spike, 484, 484) <- "A"
str_sub(omicron_spike, 493, 493) <- "R"
str_sub(omicron_spike, 496, 496) <- "S"
str_sub(omicron_spike, 498, 498) <- "R"
str_sub(omicron_spike, 501, 501) <- "Y"
str_sub(omicron_spike, 505, 505) <- "H"
str_sub(omicron_spike, 547, 547) <- "K"
str_sub(omicron_spike, 614, 614) <- "G"
str_sub(omicron_spike, 655, 655) <- "Y"
str_sub(omicron_spike, 679, 679) <- "K"
str_sub(omicron_spike, 681, 681) <- "H"
str_sub(omicron_spike, 764, 764) <- "K"
str_sub(omicron_spike, 796, 796) <- "Y"
str_sub(omicron_spike, 856, 856) <- "K"
str_sub(omicron_spike, 954, 954) <- "H"
str_sub(omicron_spike, 969, 969) <- "K"
str_sub(omicron_spike, 981, 981) <- "F"

str_sub(omicron_spike, 214, 214) <- "EPER"

str_sub(omicron_spike, 211, 212) <- ""
str_sub(omicron_spike, 141, 144) <- ""
str_sub(omicron_spike, 69, 70) <- ""


## ----omicro mhc1, eval = TRUE-------------------------------------------------
omicron_spike_peps <- mhcI(ag_present = "HLA-A*02:01", ag_stim = omicron_spike)


## ----venndiag, eval = TRUE, out.width='100%', dpi = 300-----------------------
omicron_peps_vector <- as.vector(unique(omicron_spike_peps$pep_stim))


ggVennDiagram(x=list(spike_peps_vector, omicron_peps_vector), 
              category.names= c("Original Spike Protein", "Omicron Spike Protein"), 
              set_size = 2,
              label="count") +
              theme(legend.position = "none")



## ---- eval = TRUE-------------------------------------------------------------
spike_peps %>%
  filter(!pep_stim %in% omicron_spike_peps$pep_stim) %>%
  filter(pep_stim %in% a2_spike_peps)


## ----msa prep-----------------------------------------------------------------
#  #* 1. load libraries *#
#  library(msa) # library to for alignment plots
#  library(Biostrings) # library for AAStringSet
#  library(here) # path to save plots
#  
#  #* 2. check s417 and 976 *#
#  # The K-> N mutation at 417 keeps this peptide from being a strong binder
#  s417 <- substr(spike, 417, 425) # "KIADYNYKL"
#  
#  # The L->F mutation at 981 alters this peptide's binding ability
#  s976 <- substr(spike, 976, 984) # "VLNDILSRL"
#  
#  #* 3. align spike and omicron_spike *#
#  seq_in <- AAStringSet(c(spike,  omicron_spike))
#  names(seq_in) <- c("spike", "omicron_spike")
#  aln <- msa(seq_in, type = "protein", order = "input")
#  
#  # pull out aligned sequences
#  spk_algn <- aln@unmasked$spike %>% as.character()
#  omc_algn <- aln@unmasked$omicron_spike %>% as.character()
#  
#  #* 4. find locations of s417 and s976*#
#  s417_st <- str_locate(spk_algn, s417)[1,1]
#  s417_ed <- str_locate(spk_algn, s417)[1,2]
#  s417 == str_sub(omc_algn, s417_st, s417_ed)
#  s417_algn <- str_sub(spk_algn, s417_st, s417_ed)
#  
#  s976_st <- str_locate(spk_algn, s976)[1,1]
#  s976_ed <- str_locate(spk_algn, s976)[1,2]
#  s976 == str_sub(spk_algn, s976_st, s976_ed)
#  s976_algn <- str_sub(spk_algn, s976_st, s976_ed)
#  
#  #* 5. pull out sub_sequence starting from s417 and ending at s976 *#
#  spk_algn <- str_sub(spk_algn, s417_st,  s976_ed)
#  omc_algn <- str_sub(omc_algn, s417_st,  s976_ed)
#  
#  #* 6. add - to s417 and s976 to make them have same length with aligned sub spike and sub omicron-spike *#
#  s417_n <- nchar(omc_algn) - nchar(s417_algn)
#  s417_algn <- paste0(s417_algn, strrep("-", s417_n))
#  
#  s976_n <- nchar(omc_algn) - nchar(s976_algn)
#  s976_algn <- paste0(strrep("-", s976_n), s976_algn)
#  
#  #* 6. put all sub-sequences into one AAStringSet and plot *#
#  seq4plot <- AAStringSet(c(spk_algn, omc_algn, s417_algn, s976_algn))
#  names(seq4plot) <- c("spik", "omicron_spike", "s417", "s976")
#  
#  # fake sub alignment
#  aln4plot <- as(AAMultipleAlignment(seq4plot),
#                 "MsaAAMultipleAlignment")
#  
#  fn <- "s417_"
#  seq_pos <- c(1,9)
#  
#  fn <- "s976_"
#  seq_pos <- c(nchar(s976_algn)-9+1,nchar(s976_algn))
#  
#  fn <- "s417_s976_"
#  seq_pos <- c(1,nchar(s976_algn))
#  
#  ## create PDF file according to some custom settings
#  
#  PlotAlgn <- function(fl_nm, seq_nm, pos) {
#    tmpFile <- tempfile(pattern=fl_nm,
#                        tmpdir=here(),
#                        fileext=".pdf")
#  
#    msaPrettyPrint(seq_nm,
#                   pos,
#                   file = tmpFile,
#                   output = "pdf",
#                   showNames = "left",
#                   showNumbering = "right",
#                   showLogo = "top",
#                   showConsensus = "bottom",
#                   logoColors = "rasmol",
#                   verbose = FALSE,
#                   askForOverwrite = FALSE,
#                   showLegend = FALSE)
#  }
#  
#  PlotAlgn(fl_nm = "s417_",
#           seq_nm = aln4plot,
#           pos = c(1,9))
#  
#  PlotAlgn(fl_nm = "s976_",
#           seq_nm = aln4plot,
#           pos = c(nchar(s976_algn)-9+1, nchar(s976_algn)))
#  
#  PlotAlgn(fl_nm = "s417_s976_",
#           seq_nm = aln4plot,
#           pos = c(1, nchar(s976_algn)))

## ----msa plot, eval = TRUE, out.width='100%', dpi = 300-----------------------
knitr::include_graphics(system.file("extdata/vgn/plot/covid_s417.pdf",
                                    package = "epitopeR"))

knitr::include_graphics(system.file("extdata/vgn/plot/covid_s976.pdf",
                        package = "epitopeR"))

