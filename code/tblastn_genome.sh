#!/bin/bash

#MSH6 evolution project

GC=$1
Q=$2

# BLAST databases for two species
DB=data/genomes/${GC}/genome.fna

# Query sequences for species 1
QUERY1=github/files/blast_queries/$Q

echo $QUERY1

BASENAMEQ=$(basename "$Q" .faa)

cd /home/gmonroe/projects/tol_reader_repair
mkdir -p results/tblastn/$GC
mkdir -p results/tblastn/$GC/$BASENAMEQ

# BLAST protein A from species 1 against species 2 proteome
tblastn -query $QUERY1 -db $DB -outfmt 6 > results/tblastn/$GC/$BASENAMEQ/tblastn_results.txt
