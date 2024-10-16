
#!/bin/bash -l
#MSH6 evolution project

# annotated MSH6 proteins based on interproscan and blast results
# reads in "results/predicted_msh6/all.faa.tsv",
#           "github/files/ncbi_dataset.tsv",
#           "results/predicted_msh6_blast.txt",
#            "done.txt"
# writes "results/predicted_msh6_blast_annotated.txt"
#       "results/tblastn_msh6_blast_dist.txt"
#       "results/msh6_organism_table.txt"
Rscript ./github/code/annotate_MSH6_readers.R
