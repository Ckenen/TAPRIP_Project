#!/bin/sh

# Outputs:
#   ${prefix}.bam
#   ${prefix}.bam.bai
#   ${prefix}.unmapped.fastq.1.gz
#   ${prefix}.unmapped.fastq.2.gz

fq1=$1
fq2=$2
index=$3
threads=$4
prefix=$5

bowtie2 -p ${threads} --local --no-unal -t --un-conc-gz ${prefix}.unmapped.fastq.gz -x ${index}/ref -1 ${fq1} -2 ${fq2} > ${prefix}.unsorted.sam
samtools view -@ ${threads} -u -o ${prefix}.unsorted.bam ${prefix}.unsorted.sam
samtools sort -@ ${threads} -T ${prefix}_TMP -o ${prefix}.bam ${prefix}.unsorted.bam
samtools index -@ ${threads} ${prefix}.bam
samtools flagstat -@ ${threads} > ${prefix}.flagstat
rm ${prefix}.unsorted.sam ${prefix}.unsorted.bam
