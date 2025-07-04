#!/usr/bin/env python
bam_ip=$1
bam_in=$2
bed_ip=$3
outdir=$4

mkdir -p ${outdir}
samtools view -cF 4 ${bam_ip} > ${outdir}/ip_mapped_readnum.txt
samtools view -cF 4 ${bam_in} > ${outdir}/input_mapped_readnum.txt
gzip -d -c ${bed_ip} > ${outdir}/peakClusters.bed

BINROOT="/datg/chenzonggui/paint_project/software/eclip-master/bin"

perl ${BINROOT}/overlap_peakfi_with_bam.pl \
    ${bam_ip} ${bam_in} ${outdir}/peakClusters.bed \
    ${outdir}/ip_mapped_readnum.txt \
    ${outdir}/input_mapped_readnum.txt \
    ${outdir}/peakClusters.normed.bed

perl ${BINROOT}/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.pl \
    ${outdir}/peakClusters.normed.bed \
    ${outdir}/peakClusters.normed.compressed.bed
