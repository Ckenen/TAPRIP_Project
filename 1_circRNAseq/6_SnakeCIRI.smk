#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
#SAMPLES = ["20210720_circRNAseq_oeHNRNPK_Rep1"]
SAMPLES2 = list(filter(lambda s: "_RNAseq_" in s, SAMPLES))
GENOME_FASTA = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa"
ANNOTATION_GTF = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.sorted.gtf"
FQDIR = "results/01_prepare/03_bowtie2_mapping_ribo"
OUTDIR = "results/06_ciri_quant"
CIRI_QUANT_CONFIG = "data/ciri_quant_config.yaml"

rule all:
    input:
        OUTDIR + "/01_index/bwa_index",
        OUTDIR + "/01_index/hisat2_index",
        expand(OUTDIR + "/02_ciri_quant/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/03_ciri_quant_rnaser/{sample}", sample=SAMPLES2),

rule bwa_index:
    input:
        fa = GENOME_FASTA
    output:
        directory(OUTDIR + "/01_index/bwa_index")
    log:
        OUTDIR + "/01_index/bwa_index.log"
    conda:
        "ciri"
    shell:
        """
        mkdir -p {output}
        bwa index -a bwtsw -p {output}/ref {input.fa} &> {log}
        """

rule hisat2_index:
    input:
        fa = GENOME_FASTA
    output:
        directory(OUTDIR + "/01_index/hisat2_index")
    log:
        OUTDIR + "/01_index/hisat2_index.log"
    threads:
        24
    conda:
        "ciri"
    shell:
        """
        mkdir {output}
        hisat2-build -p {threads} {input.fa} {output}/ref &> {log}
        """

rule ciri_quant:
    input:
        fq1 = FQDIR + "/{sample}.unmapped.fastq.1.gz",
        fq2 = FQDIR + "/{sample}.unmapped.fastq.2.gz",
        bwa_idx = rules.bwa_index.output,
        hs_idx = rules.hisat2_index.output,
    output:
        directory(OUTDIR + "/02_ciri_quant/{sample}")
    log:
        OUTDIR + "/02_ciri_quant/{sample}.log"
    threads:
        24
    conda:
        "CIRIquant_env"
    shell:
        """
        CIRIquant -t 12 -l 2 -o {output} -p {wildcards.sample} \
            --config {CIRI_QUANT_CONFIG} -1 {input.fq1} -2 {input.fq2} &> {log}
        """

def get_rnaser_sample(sample):
    assert "_RNAseq_" in sample
    return sample.replace("_RNAseq_", "_circRNAseq_")

rule ciri_quant_rnaser:
    input:
        fq1 = FQDIR + "/{sample}.unmapped.fastq.1.gz",
        fq2 = FQDIR + "/{sample}.unmapped.fastq.2.gz",
        bwa_idx = rules.bwa_index.output,
        hs_idx = rules.hisat2_index.output,
        rnaser = lambda wildcards: OUTDIR + "/02_ciri_quant/%s" % get_rnaser_sample(wildcards.sample)
    output:
        directory(OUTDIR + "/03_ciri_quant_rnaser/{sample}")
    log:
        OUTDIR + "/03_ciri_quant_rnaser/{sample}.log"
    params:
        rnaser_sample = lambda wildcards: get_rnaser_sample(wildcards.sample)
    threads:
        12
    conda:
        "CIRIquant_env"
    shell:
        """
        CIRIquant -t 12 -l 2 -o {output} -p {wildcards.sample} \
            --RNaseR {input.rnaser}/{params.rnaser_sample}.gtf \
            --config {CIRI_QUANT_CONFIG} -1 {input.fq1} -2 {input.fq2} &> {log}
        """