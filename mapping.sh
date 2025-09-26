# !/bin/bash

# Rename to A_thaliana.fa
mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa A_thaliana.fa

# create genome index directory
STAR --runMode genomeGenerate --genomeDir genomeIndex --genomeFastaFiles A_thaliana.fa

# create directories for output
mkdir -p mapped

for infile in trimmed/*.fastq.gz ; do

	gunzip $infile

done

for infile in trimmed/*.fastq ; do
	outfile=$(basename $infile .fastq)
	STAR --genomeDir genome/genomeIndex --readFilesIn $infile --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done

for file in mapped/*.bam; do

   samtools view $file | head -n 3
   samtools view -c $file | head -n 3
   samtools  flagstat $file | head -n 3

done

# create a directory for IGV viewing
mkdir IGV
cp mapped/*.bam IGV/
