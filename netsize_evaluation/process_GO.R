###########################################################
# This script downloads the GO terms from pseudomonas.com and processes it
# into the right format.
#
# Usage:
#     Rscript process_GO.R out_file 
#
#     out_file: file path to save the GO term output
#
###########################################################

pacman::p_load("readr", "dplyr")

out_file <- commandArgs(trailingOnly = TRUE)[1]

GO_url <- "http://pseudomonas.com/goterms/list?goEvidenceCode=MC&format=TAB&extension=.txt&strain_id=107&queryStrain=Pseudomonas%20aeruginosa%20PAO1%20(Reference)"
# directly read in the data from pseudomonas.com
GO <- read_tsv(GO_url)
# only consider biological processes
BP <- GO %>% filter(Namespace == "biological_process")
# built a GO term ID column
BP <- BP %>% mutate(GO_name = paste(Accession, `GO Term`, sep = "_"))
# collapse genes by GO term
BP_summary <- BP %>% group_by(GO_name) %>%
  summarise(gene_count = n(), genes = paste(`Locus Tag`, collapse = ";"))
# only consider pathways with more than 5 genes and less than 100 genes as
# meaningful pathways
BP_summary <- BP_summary %>% filter(gene_count >= 5 & gene_count <= 100)

readr::write_tsv(BP_summary, out_file, col_names = FALSE)