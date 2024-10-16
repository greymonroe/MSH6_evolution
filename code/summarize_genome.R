# Load required libraries
library(data.table)
library(seqinr)

# Source the required functions
source("code/summarize_gff.R")
source("code/summarize_amino_acids.R")

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) != 4) {
  stop("Please provide four arguments: gff_file, gff_out_file, protein_file, protein_out_file")
}

# Assign arguments to variables
gff_file <- args[1]
gff_out_file <- args[2]
protein_file <- args[3]
protein_out_file <- args[4]

# Call the functions
summarize_gff(gff_file, gff_out_file)
summarize_amino_acids(protein_file, protein_out_file)
