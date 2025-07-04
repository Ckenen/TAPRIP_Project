#!/bin/sh
set -e

bam_rep1_ip=$1
bam_rep2_ip=$2
bam_rep1_in=$3
bam_rep2_in=$4
bed_rep1_ip=$5
bed_rep2_ip=$6
outdir=$7

mkdir -p ${outdir}
samtools view -c -F 4 ${bam_rep1_ip} > ${outdir}/rep1_clip.readnum
samtools view -c -F 4 ${bam_rep2_ip} > ${outdir}/rep2_clip.readnum
samtools view -c -F 4 ${bam_rep1_in} > ${outdir}/rep1_input.readnum
samtools view -c -F 4 ${bam_rep2_in} > ${outdir}/rep2_input.readnum
gzip -d -c ${bed_rep1_ip} > ${outdir}/rep1_peakClusters.bed
gzip -d -c ${bed_rep2_ip} > ${outdir}/rep2_peakClusters.bed

BINROOT="/datg/chenzonggui/paint_project/software/merge_peaks-master/bin"

perl ${BINROOT}/perl/overlap_peakfi_with_bam.pl \
    ${bam_rep1_ip} \
    ${bam_rep1_in} \
    ${outdir}/rep1_peakClusters.bed \
    ${outdir}/rep1_clip.readnum \
    ${outdir}/rep1_input.readnum \
    ${outdir}/rep1_normed_peaks.bed

perl ${BINROOT}/perl/overlap_peakfi_with_bam.pl \
    ${bam_rep2_ip} \
    ${bam_rep2_in} \
    ${outdir}/rep2_peakClusters.bed \
    ${outdir}/rep2_clip.readnum \
    ${outdir}/rep2_input.readnum \
    ${outdir}/rep2_normed_peaks.bed

perl ${BINROOT}/perl/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl \
    ${outdir}/rep1_normed_peaks.bed.full \
    ${outdir}/rep1_normed_peaks.compressed.bed \
    ${outdir}/rep1_normed_peaks.compressed.bed.full

perl ${BINROOT}/perl/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl \
    ${outdir}/rep2_normed_peaks.bed.full \
    ${outdir}/rep2_normed_peaks.compressed.bed \
    ${outdir}/rep2_normed_peaks.compressed.bed.full

perl ${BINROOT}/perl/make_informationcontent_from_peaks.pl \
    ${outdir}/rep1_normed_peaks.compressed.bed.full \
    ${outdir}/rep1_clip.readnum \
    ${outdir}/rep1_input.readnum \
    ${outdir}/rep1_normed_peaks.compressed.bed.entropy.full \
    ${outdir}/rep1_normed_peaks.compressed.bed.entropy.excessreads

perl ${BINROOT}/perl/make_informationcontent_from_peaks.pl \
    ${outdir}/rep2_normed_peaks.compressed.bed.full \
    ${outdir}/rep2_clip.readnum \
    ${outdir}/rep2_input.readnum \
    ${outdir}/rep2_normed_peaks.compressed.bed.entropy.full \
    ${outdir}/rep2_normed_peaks.compressed.bed.entropy.excessreads

python ${BINROOT}/full_to_bed.py \
    --input ${outdir}/rep1_normed_peaks.compressed.bed.entropy.full \
    --output ${outdir}/rep1_normed_peaks.compressed.bed.entropy.bed

python ${BINROOT}/full_to_bed.py \
    --input ${outdir}/rep2_normed_peaks.compressed.bed.entropy.full \
    --output ${outdir}/rep2_normed_peaks.compressed.bed.entropy.bed

source /home/chenzonggui/miniconda3/etc/profile.d/conda.sh
conda activate idr
idr --samples \
    ${outdir}/rep1_normed_peaks.compressed.bed.entropy.bed \
    ${outdir}/rep2_normed_peaks.compressed.bed.entropy.bed \
    --input-file-type bed \
    --rank 5 \
    --peak-merge-method max \
    --plot \
    -o ${outdir}/01v02.idr.out
conda deactivate

perl ${BINROOT}/perl/parse_idr_peaks.pl \
    ${outdir}/01v02.idr.out \
    ${outdir}/rep1_normed_peaks.compressed.bed.entropy.full \
    ${outdir}/rep2_normed_peaks.compressed.bed.entropy.full \
    ${outdir}/01v02.idr.out.bed
