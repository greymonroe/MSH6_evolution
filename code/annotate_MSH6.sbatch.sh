#!/bin/bash -l
#SBATCH -o /home/gmonroe/slurm-log2/%j-stdout.txt
#SBATCH -e /home/gmonroe/slurm-log2/%j-stderr.txt  # Redirecting stderr to a separate file
#SBATCH -J annotate
#SBATCH -t 1-00:00:00
#SBATCH --partition=bmh
#SBATCH --ntasks=32  # Number of tasks (cores) set to 16
#SBATCH --mem=128G    # Memory set to 32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=greymonroe@gmail.com

#!/bin/bash -l

# annotated MSH6 proteins based on interproscan and blast results
# reads in "results/predicted_msh6/all.faa.tsv", 
#           "github/files/ncbi_dataset.tsv", 
#           "results/predicted_msh6_blast.txt",
#            "done.txt"
# writes "results/predicted_msh6_blast_annotated.txt"
#       "results/tblastn_msh6_blast_dist.txt"
#       "results/msh6_organism_table.txt"
Rscript ./github/code/annotate_MSH6_readers.R