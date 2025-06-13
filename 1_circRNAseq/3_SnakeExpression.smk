#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
BAMDIR = "results/02_mapping/04_rmdup"
OUTDIR = "results/03_expression"
#samples = samples[:-16]

rule all:
    input:
        expand(OUTDIR + "/01_stringtie2/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/02_fpkm/{sample}.tsv", sample=SAMPLES),

rule stringtie2:
    input:
        bam = BAMDIR + "/{sample}.bam",
        gtf = config["GENOME_GTF"]
    output:
        out = directory(OUTDIR + "/01_stringtie2/{sample}")
    log:
        OUTDIR + "/01_stringtie2/{sample}.log"
    threads:
        8
    shell:
        """
        mkdir {output}
        stringtie {input.bam} --rf -e -B -p {threads} -A {output.out}/gene_abund.tab -C {output.out}/cov_refs.gtf \
            -G {input.gtf} -o {output.out}/transcripts.gtf &> {log}
        pigz -p {threads} {output.out}/*
        """

rule calculate_fpkm:
    input:
        bam = BAMDIR + "/{sample}.bam",
        bed = config["TRANSCRIPT_BED_GZ"],
        tsv = config["ANNOTATION_TSV"]
    output:
        tsv = OUTDIR + "/02_fpkm/{sample}.tsv"
    log:
        OUTDIR + "/02_fpkm/{sample}.log"
    threads:
        12
    shell:
        """
        ./scripts/expression/calculate_fpkm.py {input.bam} {input.bed} {input.tsv} {threads} {output.tsv}
        """

