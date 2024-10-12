#!/bin/bash -l 

################################################################################
# Script Name: download_ncbi.sh
# Description: This script downloads specific file types for a given genome 
#              accession (GC value) from the NCBI database. After downloading, 
#              BLAST databases are generated for the protein and genomic sequences.
#
# File Types Downloaded from NCBI:
#   - Protein sequences (PROT_FASTA)
#   - Genomic sequences (GENOME_FASTA)
#   - Genome annotations (GENOME_GFF)
#
# BLAST Databases Generated:
#   - Protein sequences (protein.faa)
#   - Genomic sequences (genome.fna)
#
# Directory Structure:
#   Input:
#     - None (Files are downloaded from NCBI)
#
#   Output:
#     - Genomic and protein sequences, annotations, and BLAST databases are 
#       stored in: data/genomes/[GC_VALUE]/
#
# Usage:
#   ./download_ncbi.sh [GC_VALUE]
#
# Parameters:
#   GC_VALUE - The GC value from the NCBI database identifying the genome.
#
# Example:
#   ./download_ncbi.sh GC12345
#
# Notes:
#   - The script uses the curl tool to download files from NCBI.
#   - If the download or unzip operations fail, the script will retry until successful.
#   - The script activates the 'tol_reader_repair' conda environment for operations.
################################################################################


conda activate tol_reader_repair

GC=$1

# Create a directory with the name of the accession
mkdir -p data/genomes/$GC
mkdir -p results/genome_stats/$GC

# Set the output directory variable
outdir=data/genomes/$GC

# Navigate to the output directory
cd $outdir

rm -rf *

SUCCESS=false

# Keep trying until we've downloaded a successful zip file
until $SUCCESS; do

    
while ! curl --speed-time 30 --speed-limit 12500 -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/${GC}/download?include_annotation_type=PROT_FASTA&include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&filename=${GC}.zip" -H "Accept: application/zip"; do
    echo "Retrying download..."
    rm *
    sleep 5  # wait for 5 seconds before retrying
done
# Download the file using curl to the current directory
    #curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/${GC}/download?include_annotation_type=PROT_FASTA&include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&filename=${GC}.zip" -H "Accept: application/zip"

    # Attempt to unzip the file
    if unzip ${GC}.zip; then
        SUCCESS=true
    else
        echo "Unzip failed. Retrying download..."
        rm -rf *
    fi

done

# Move the desired files to the top-level directory
mv ncbi_dataset/data/$GC/* .

# Check if protein.faa is present before running makeblastdb
if [ -f protein.faa ]; then
    makeblastdb -in protein.faa -dbtype prot
else
    echo "protein.faa not found!"
fi

# Check if there's exactly one file with the pattern *genomic.fna before copying
if [ $(ls *genomic.fna 2>/dev/null | wc -l) -eq 1 ]; then
    mv *genomic.fna genome.fna
else
    echo "Either no *genomic.fna file found or multiple files found!"
fi

# Check if genome.fna is present before running makeblastdb
if [ -f genome.fna ]; then
    makeblastdb -in genome.fna -dbtype nucl
else
    echo "genome.fna not found!"
fi


# Optional: If you want to clean up, you can remove the ncbi_dataset directory
rm -rf ncbi_dataset/

if [ -f genomic.gff ] && [ -f protein.faa ]; then
echo getting genome stats
Rscript ~/projects/tol_reader_repair/github/code/summarize_genome.R genomic.gff ~/projects/tol_reader_repair/results/genome_stats/$GC/CDS_stats.csv protein.faa ~/projects/tol_reader_repair/results/genome_stats/$GC/AA_freq.csv
fi

cd /home/gmonroe/projects/tol_reader_repair

conda deactivate