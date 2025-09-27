#!/bin/bash
set -euo pipefail  # exit on errors, undefined vars, and broken pipes
trap 'echo "Error at line $LINENO"; exit 1' ERR

# Check if input reference exists
if [[ ! -f Arabidopsis_thaliana.TAIR10.dna.toplevel.fa ]]; then
    echo "Reference FASTA not found."
    exit 1
fi

# Rename to A_thaliana.fa (only if not already renamed)
if [[ ! -f A_thaliana.fa ]]; then
    mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa A_thaliana.fa
fi

# Create genome index directory
mkdir -p genomeIndex
if ! STAR --runMode genomeGenerate --genomeDir genomeIndex --genomeFastaFiles A_thaliana.fa; then
    echo "STAR genomeGenerate failed."
    exit 1
fi

# Create directories for output
mkdir -p mapped

# Decompress fastq.gz safely
for infile in trimmed/*.fastq.gz ; do
    if [[ -f "$infile" ]]; then
        gunzip -k "$infile"  # -k keeps original gz file
    else
        echo "No .fastq.gz files found in trimmed/"
        exit 1
    fi
done

# Run STAR mapping
for infile in trimmed/*.fastq ; do
    outfile=$(basename "$infile" .fastq)
    if ! STAR --genomeDir genomeIndex --readFilesIn "$infile" \
        --outFileNamePrefix mapped/"$outfile" \
        --outSAMtype BAM SortedByCoordinate --outSAMattributes All; then
        echo "STAR mapping failed for $infile"
        exit 1
    fi
done

# Samtools checks
for file in mapped/*.bam; do
    if [[ -f "$file" ]]; then
        samtools view "$file" | head -n 3 || echo "samtools view failed on $file"
        samtools view -c "$file" || echo "samtools count failed on $file"
        samtools flagstat "$file" | head -n 3 || echo "samtools flagstat failed on $file"
    else
        echo "No BAM files found in mapped/"
        exit 1
    fi
done

# Prepare IGV directory
mkdir -p IGV
cp mapped/*.bam IGV/ || echo "No BAM files to copy into IGV/"
