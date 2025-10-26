#!/bin/sh

# STAR index for human (v2.7.9a)
mkdir -p star_hg38_index
STAR --runMode genomeGenerate --runThreadN 12 --genomeDir star_hg38_index --genomeFastaFiles ucsc_hg38/hg38.fa --sjdbGTFfile ucsc_hg38/knownGene.gtf &> star_hg38_index.log

# STAR index for human repbase (v2.7.9a)
mkdir -p star_humrep
STAR --runMode genomeGenerate --runThreadN 4 --genomeDir star_humrep --genomeFastaFiles RepBase/RepBase26.08.fasta/humrep.ref &> star_humrep.log

