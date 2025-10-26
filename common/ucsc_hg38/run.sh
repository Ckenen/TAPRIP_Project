#!/bin/sh

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -d -c hg38.fa.gz > hg38.fa
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz
zcat knownGene.txt.gz | awk -v FS='\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | genePredToGtf file stdin stdout | bedtools sort -i - > knownGene.gtf
zcat knownGene.txt.gz | awk -v FS='\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | genePredToBed stdin stdout | bedtools sort -i - > knownGene.bed
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
gzip -d -c refFlat.txt.gz > refFlat.txt
gzip -d -c knownGene.txt.gz > knownGene.txt