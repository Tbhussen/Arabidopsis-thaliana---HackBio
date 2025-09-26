#!/bin/bash

# Download annotation file
#wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz

#unzip the gff3 file
#gunzip Arabidopsis_thaliana.TAIR10.62.gff3.gz

#rename the gff3 file to c_elegans.gff3
# this is the file that will be used for featureCounts
#mv Arabidopsis_thaliana.TAIR10.62.gff3 A_thaliana.gff3

#run featureCounts
# -O: count reads overlapping multiple features
# -t: feature type to count (gene)
# -g: attribute to group features (ID)
# -a: annotation file (GFF3)
# -o: output file

#mkdir -p counts

echo Apply Feature Counts to the BAM files

featureCounts -O -t gene -g ID -a genome/A_thaliana.gff3 -o counts/counts.txt IGV/*.bam

cd counts/
