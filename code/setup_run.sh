#!/bin/bash -l

conda activate tol_reader_repair

mkdir data
mkdir data/genomes
mkdir data/references
mkdir results
mkdir results/blastp
mkdir results/blastn
mkdir results/blastx
mkdir results/tblastn
mkdir results/genome_stats
mkdir ncbi_batches

./github/code/download_ncbi.sh GCA_000001735.2
mv data/genomes/GCA_000001735.2 data/references/
  seqkit seq -n data/references/GCA_000001735.2/protein.faa > data/references/GCA_000001735.2/protein_names.txt
./github/code/download_ncbi.sh GCF_000001405.40
mv data/genomes/GCF_000001405.40 data/references/
  seqkit seq -n data/references/GCF_000001405.40/protein.faa > data/references/GCF_000001405.40/protein_names.txt
conda deactivate

# Large batch command (96 jobs) to run in parallel across genomes
### Identifying MSH6 orthologs and potential reader domains ###
# downloads NCBI genomes and proteomes for all species and run blastp and tblastn
./github/code/batches.sh 96 ./github/code/find_ortho_MSH6_proteins.sbatch.sh

# predict MSH6 orthologs based on reciprocal blastp results and annotate proteins with interproscan
# runs ./github/code/identify_predicted_MSH6_orthologs.R
# runs interproscan
sbatch ./github/code/interproscan_predicted_MSH6.sbatch.sh

# parse results from tblastn, blastp, interproscan to annotate predicted MSH6 orthologs
# final result are files with each protein annotated with direct evidence of Tudor or PWWP domains
# also identifies potential mis-annotations (calculates distance between tblastn hits for MSH6 and domain) which will be classified as ambiguous
sbatch ./github/code/annotate_MSH6.sbatch.sh


