#!bin/bash

mkdir -p qc_reports

for file in raw_data/*.fastq.gz; do

   fastqc $file -o qc_reports/

done

#Trimming
mkdir -p trimmed

for file in raw_data/*.fastq.gz; do

   sample=$(basename $file ".fastq.gz")
   fastp -i $file -o trimmed/${sample}.trim.fastq.gz -h trimmed/${sample}.html

done

echo The quality control is Complete
