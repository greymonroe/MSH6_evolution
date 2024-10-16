
#!/bin/bash -l

BATCH_FILE=$1
conda activate tol_reader_repair
touch batches_progress/$SLURM_JOB_ID.txt
touch batches_progress_missing_data/${SLURM_JOB_ID}.txt

while read -r GC; do
    echo $GC
    # Download
    ./github/code/download_ncbi.sh $GC

    # Check if the necessary files are present
    if [ -f data/genomes/$GC/genome.fna ] && [ -f data/genomes/$GC/protein.faa ]; then
        ./github/code/recip_blastp_protein.sh GCF_000001405.40 $GC human_msh6.faa
        ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_msh6.faa
        ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_msh6_tudor.faa
         ./github/code/recip_blastp_protein.sh GCF_000001405.40 $GC human_msh6_pwwp.faa
        ./github/code/tblastn_genome.sh $GC human_msh6.faa
         ./github/code/tblastn_genome.sh $GC arabidopsis_msh6.faa
         ./github/code/tblastn_genome.sh $GC human_msh6_pwwp.faa
         ./github/code/tblastn_genome.sh $GC arabidopsis_msh6_tudor.faa
        #tar -czvf data/genomes/$GC.tar.gz data/genomes/$GC/
        echo $GC >> batches_progress/$SLURM_JOB_ID.txt
    else
        rm -rf data/genomes/$GC
        echo $GC
        echo $GC >> batches_progress_missing_data/${SLURM_JOB_ID}.txt
    fi
done < $BATCH_FILE


