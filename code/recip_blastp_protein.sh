#!/bin/bash

GC1=$1
GC2=$2
Q=$3

# BLAST databases for two species
DB1=data/references/${GC1}/protein.faa
DB2=data/genomes/${GC2}/protein.faa

# Query sequences for species 1
QUERY1=github/files/blast_queries/$Q

echo $QUERY1

BASENAMEQ=$(basename "$Q" .faa)

cd /home/gmonroe/projects/tol_reader_repair
mkdir -p results/blastp/$GC2
mkdir -p results/blastp/$GC2/$BASENAMEQ
mkdir -p results/blastp/$GC2/$BASENAMEQ/proteins
mkdir -p results/blastp/$GC2/$BASENAMEQ/rev_blastp
mkdir -p results/blastp/$GC2/$BASENAMEQ/gff


# Create BLAST databases if they don't exist
#makeblastdb -in $DB1 -dbtype prot
#makeblastdb -in $DB2 -dbtype prot

# BLAST protein A from species 1 against species 2 proteome
blastp -query $QUERY1 -db $DB2 -outfmt 6 -evalue 0.01 > results/blastp/$GC2/$BASENAMEQ/blastp_results.txt

awk '{print $2}' results/blastp/$GC2/$BASENAMEQ/blastp_results.txt > results/blastp/$GC2/$BASENAMEQ/top_hits_ids.txt

# Loop over top hits and BLAST against species 1 proteome
while read -r TOP_HIT; do
    
    grep $TOP_HIT data/genomes/${GC2}/genomic.gff > results/blastp/$GC2/$BASENAMEQ/gff/${TOP_HIT}.gff

    # Extract the sequence
    seqkit grep -r -p "$TOP_HIT" $DB2 > results/blastp/$GC2/$BASENAMEQ/proteins/${TOP_HIT}.fasta

    # BLAST against species 1 proteome
    blastp -query results/blastp/$GC2/$BASENAMEQ/proteins/${TOP_HIT}.fasta -db $DB1 -evalue 1e-5 -outfmt 6 > results/blastp/$GC2/$BASENAMEQ/rev_blastp/${TOP_HIT}.txt
    
done < results/blastp/$GC2/$BASENAMEQ/top_hits_ids.txt

