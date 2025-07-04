#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
BAMDIR = "results/01_pipeline/09_dedup"
OUTDIR = "results/03_tracks"

rule all:
    input:
        expand(OUTDIR + "/01_bw/{sample}_raw_both.bw", sample=SAMPLES),
        expand(OUTDIR + "/02_bin_coverage/{sample}.txt", sample=SAMPLES),

rule bam2bw:
    input:
        bam = BAMDIR + "/{sample}.dedup.sorted.bam"
    output:
        bw = OUTDIR + "/01_bw/{sample}_raw_both.bw"
    log:
        OUTDIR + "/01_bw/{sample}.log"
    params:
        prefix = OUTDIR + "/01_bw/{sample}"
    shell:
        """
        ./scripts/bam2bw.sh {input.bam} {params.prefix} &> {log}
        """

rule stat_bin_coverage:
    input:
        bw = rules.bam2bw.output.bw
    output:
        txt = OUTDIR + "/02_bin_coverage/{sample}.txt"
    log:
        OUTDIR + "/02_bin_coverage/{sample}.log"
    threads:
        12
    shell:
        """
        ./scripts/stat_bin_coverage.py {input.bw} {threads} {output.txt} &> {log}
        """