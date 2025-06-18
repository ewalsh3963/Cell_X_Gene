# /usr/bin/R 
# Evan Walsh (evan.walsh@capsida.com)

# options(width = Sys.getenv("COLUMNS"))

## Install packages as needed
httr::set_config(httr::config(ssl_verifypeer = FALSE))
.cran_packages <- c("tidyverse", 'ggpubr', 'grid','gridExtra')
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) install.packages(.cran_packages[!.inst], repos='http://cran.us.r-project.org')
sapply(.cran_packages, require, character.only = TRUE)



.bioc_packages <- c('biomaRt')
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) BiocManager::install(.bioc_packages[!.inst])

## Load packages
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

## I/O
args = commandArgs(trailingOnly = T) 
indir = args[1]
study = args[2]
species = args[3]

indir <- '~/scratch/scanpy/MULTI_TEST'
study <- 'MULTI_TEST'
species <- 'Mus musculus'

## read in the csv file 
dat_file = file.path(indir, paste0(study, "_vars.csv"))
dat = read.table(dat_file, sep=',', header = F, stringsAsFactors = F)
colnames(dat) = 'Symbol'
gene_symbols = dat$Symbol


# Split the string by "_" and modify
parts <- unlist(strsplit(species, " "))  # Splits into c("Macaca", "Fals")
species_biomart <- paste0(tolower(substr(parts[1], 1, 1)), tolower(parts[2]))
ensembl=useMart("ensembl")
instance_mart <- useMart('ensembl', dataset = paste0(species_biomart, '_gene_ensembl'))

filter = "external_gene_name"
attribute = c("external_gene_name", "ensembl_gene_id", "hsapiens_homolog_associated_gene_name", "hsapiens_homolog_ensembl_gene")

## query for ensembl IDs
dat = getBM(attributes=attribute, 
      filters = filter,
      values = gene_symbols, 
      mart = instance_mart)

## save
basename = paste0(study, "_homo_sapiens_vars.csv")
outfile = file.path(indir, basename)
print(outfile)
write.table(dat, file = outfile, sep = ',', col.names = T, row.names = F, quote = F)