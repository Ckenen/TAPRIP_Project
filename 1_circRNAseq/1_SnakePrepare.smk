#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
FQDIR = "data/datasets"
OUTDIR = "results/1_prepare"

rule all:
    input:
        expand(OUTDIR + "/1_cutadapt/{sample}_R1.fastq.gz", sample=SAMPLES),
        OUTDIR + "/2_index/ribo_index",
        expand(OUTDIR + "/3_bowtie2/{sample}.bam", sample=SAMPLES),

# Preprocessing

rule cutadapt:
    input:
        fq1 = FQDIR + "/{sample}_R1.fastq.gz",
        fq2 = FQDIR + "/{sample}_R2.fastq.gz"
    output:
        fq1 = OUTDIR + "/1_cutadapt/{sample}_R1.fastq.gz",
        fq2 = OUTDIR + "/1_cutadapt/{sample}_R2.fastq.gz"
    log:
        OUTDIR + "/1_cutadapt/{sample}.log"
    threads:
        12
    conda:
        "cutadapt"
    shell:
        """
        cutadapt -j {threads} -m 20 -q 30 -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
            -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} &> {log}
        """

rule bowtie2_build:
    input:
        fa = config["RIBOSOME"]
    output:
        out = directory(OUTDIR + "/2_index/ribo_index")
    log:
        OUTDIR + "/2_index/ribo_index.log"
    threads:
        4
    conda:
        "bowtie2"
    shell:
        """
        mkdir -p {output.out}
        bowtie2-build --threads {threads} {input.fa} {output.out}/ref &> {log}
        """

rule bowtie2:
    input:
        fq1 = rules.cutadapt.output.fq1,
        fq2 = rules.cutadapt.output.fq2,
        idx = rules.bowtie2_build.output.out,
    output:
        bam = OUTDIR + "/3_bowtie2/{sample}.bam",
        bai = OUTDIR + "/3_bowtie2/{sample}.bam.bai",
        fq1 = OUTDIR + "/3_bowtie2/{sample}.unmapped.fastq.1.gz",
        fq2 = OUTDIR + "/3_bowtie2/{sample}.unmapped.fastq.2.gz",
    log:
        OUTDIR + "/3_bowtie2/{sample}.log"
    params:
        prefix = OUTDIR + "/3_bowtie2/{sample}"
    threads:
        12
    conda:
        "bowtie2"
    shell:
        """
        ./scripts/bowtie2_keep_unmapped_fastq_pe.sh {input} {threads} {params.prefix} &> {log}
        """
