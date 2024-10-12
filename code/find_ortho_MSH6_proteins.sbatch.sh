#!/bin/bash -l
#SBATCH -o /home/gmonroe/slurm-log2/%j-stdout.txt
#SBATCH -e /home/gmonroe/slurm-log2/%j-stderr.txt  # Redirecting stderr to a separate file
#SBATCH -J blast
#SBATCH -t 14-00:00:00
#SBATCH --partition=bmh
#SBATCH --ntasks=1  
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=greymonroe@gmail.com

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
        #./github/code/recip_blastp_protein.sh GCF_000001405.40 $GC human_msh6.faa
       # ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_msh6.faa
        #./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_msh6_tudor.faa
         ./github/code/recip_blastp_protein.sh GCF_000001405.40 $GC human_msh6_pwwp.faa
    #    ./github/code/tblastn_genome.sh $GC human_msh6.faa
    #     ./github/code/tblastn_genome.sh $GC arabidopsis_msh6.faa
         ./github/code/tblastn_genome.sh $GC human_msh6_pwwp.faa
    #     ./github/code/tblastn_genome.sh $GC arabidopsis_msh6_tudor.faa
    # ./github/code/recip_blastp_protein.sh GCF_000001405.40 $GC human_ledgfp75.faa
    # ./github/code/recip_blastp_protein.sh GCF_000001405.40 $GC human_pds5a.faa
    #     ./github/code/recip_blastp_protein.sh GCF_000001405.40 $GC human_pds5b.faa
    #     ./github/code/recip_blastp_protein.sh GCF_000001405.40 $GC human_ZCWPW1.faa
    #     ./github/code/recip_blastp_protein.sh GCF_000001405.40 $GC human_ZCWPW2.faa
    # ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_pds5a.faa
    # ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_pds5b.faa
    # ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_pds5c.faa
    # ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_pds5d.faa
    # ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_pds5e.faa
    #  # ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_atx1_tudor.faa
    #  # ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_atx2_tudor.faa
    #    ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_pds5a_tudor.faa
    #     ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_pds5b_tudor.faa
    #      ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_pds5c_tudor.faa
    #      ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_pds5d_tudor.faa
    #      ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_pds5e_tudor.faa
      #  ./github/code/recip_blastp_protein.sh GCA_000001735.2 $GC arabidopsis_AT5G10950_tudor.faa
        rm -rf data/genomes/$GC
        #tar -czvf data/genomes/$GC.tar.gz data/genomes/$GC/
        echo $GC >> batches_progress/$SLURM_JOB_ID.txt
    else
        rm -rf data/genomes/$GC
        echo $GC
        echo $GC >> batches_progress_missing_data/${SLURM_JOB_ID}.txt
    fi
done < $BATCH_FILE


