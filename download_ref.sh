#!/bin/bash

# Set variables
GENOME_URL="ftp://ftp.ensemblgenomes.org/pub/plants/release-61/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
GENOME_FILE="Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"

# Download the genome
echo "Downloading Arabidopsis thaliana genome..."
wget $GENOME_URL

# Unzip the genome file
echo "Unzipping genome file..."
gunzip $GENOME_FILE

# The FASTA file after unzipping
FASTA_FILE="Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"

# Show the first 10 lines of the FASTA file
echo "First 10 lines of the genome FASTA:"
head $FASTA_FILE

# Count the number of chromosomes/scaffolds
echo "Number of sequences (chromosomes/scaffolds) in the FASTA:"
grep -c '^>' $FASTA_FILE
