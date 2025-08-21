# Load necessary libraries if not already installed

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressMessages(library("optparse"))

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to disease GWAS summary file (SNP, A1, A2, Z) [required]"),
  make_option("--pos_dir", action="store", default=NA, type='character',
              help="*.pos files folder [required]"),  
  make_option("--model_dir", action="store", default=NA, type='character',
              help="*.RDat files folder [required]"),
  make_option("--ref_ld_dir", action="store", default=NA, type='character',
              help="pre extracted LD reference data [required]"),
  make_option("--TSS", action="store", default=NA, type='character',
              help="TSS file [required]"),                
  make_option("--output_dir", action="store", default=NA, type='character',
              help="Path to output files [required]")

)

# Get command-line arguments
opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$sumstats) || is.na(opt$pos_dir) || is.na(opt$model_dir) || is.na(opt$ref_ld_dir)|| is.na(opt$output_dir) ) {
    cat("Warning: One or more required options are missing.\n")
    q()
}

source("00_analysis_code/support.R")





########################
### association test ###
########################
association_test(
    sumstats = opt$sumstats, 
    pos_dir = opt$pos_dir, 
    model_dir = opt$model_dir, 
    ref_ld_dir = opt$ref_ld_dir, 
    TSS_file = opt$TSS,
    output_dir = opt$output_dir)
