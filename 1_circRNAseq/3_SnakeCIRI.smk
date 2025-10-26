#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
SAMPLES2 = list(filter(lambda s: "_RNAseq_" in s, SAMPLES))
CIRI_QUANT_CONFIG = "data/ciri_quant_config.yaml"
OUTDIR = "results/3_ciri"

rule all:
    input:
        OUTDIR + "/1_index/bwa_index",
        OUTDIR + "/1_index/hisat2_index",
        expand(OUTDIR + "/2_ciri_quant/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/3_ciri_quant_rnaser/{sample}", sample=SAMPLES2),

rule bwa_index:
    input:
        fa = "/home/zgchen/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa"
    output:
        out = directory(OUTDIR + "/1_index/bwa_index")
    log:
        OUTDIR + "/1_index/bwa_index.log"
    conda:
        "ciri"
    shell:
        """
        mkdir -p {output}
        bwa index -a bwtsw -p {output}/ref {input.fa} &> {log}
        """

rule hisat2_index:
    input:
        fa = "/home/zgchen/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa"
    output:
        out = directory(OUTDIR + "/1_index/hisat2_index")
    log:
        OUTDIR + "/1_index/hisat2_index.log"
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
        fq1 = "results/1_prepare/2_bowtie2/{sample}.unmapped.fastq.1.gz",
        fq2 = "results/1_prepare/2_bowtie2/{sample}.unmapped.fastq.2.gz",
        bwa_idx = rules.bwa_index.output,
        hs_idx = rules.hisat2_index.output,
    output:
        out = directory(OUTDIR + "/2_ciri_quant/{sample}")
    log:
        OUTDIR + "/2_ciri_quant/{sample}.log"
    threads:
        24
    conda:
        "CIRIquant_env"
    shell:
        """
        CIRIquant -t {threads} -l 2 -o {output} -p {wildcards.sample} \
            --config {CIRI_QUANT_CONFIG} -1 {input.fq1} -2 {input.fq2} &> {log}
        """

def get_rnaser_sample(sample):
    assert "_RNAseq_" in sample
    return sample.replace("_RNAseq_", "_circRNAseq_")

rule ciri_quant_rnaser:
    input:
        fq1 = "results/1_prepare/2_bowtie2/{sample}.unmapped.fastq.1.gz",
        fq2 = "results/1_prepare/2_bowtie2/{sample}.unmapped.fastq.2.gz",
        bwa_idx = rules.bwa_index.output,
        hs_idx = rules.hisat2_index.output,
        rnaser = lambda wildcards: rules.ciri_quant.output.out.format(sample=get_rnaser_sample(wildcards.sample))
    output:
        out = directory(OUTDIR + "/3_ciri_quant_rnaser/{sample}")
    log:
        OUTDIR + "/3_ciri_quant_rnaser/{sample}.log"
    params:
        rnaser_sample = lambda wildcards: get_rnaser_sample(wildcards.sample)
    threads:
        12
    conda:
        "CIRIquant_env"
    shell:
        """
        CIRIquant -t {threads} -l 2 -o {output} -p {wildcards.sample} \
            --RNaseR {input.rnaser}/{params.rnaser_sample}.gtf \
            --config {CIRI_QUANT_CONFIG} -1 {input.fq1} -2 {input.fq2} &> {log}
        """
