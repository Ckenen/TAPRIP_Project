#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
BAMDIR = "results/1_pipeline/9_dedup"
OUTDIR = "results/3_tracks"

rule all:
    input:
        expand(OUTDIR + "/1_bw/{sample}_raw_both.bw", sample=SAMPLES),

rule bam2bw:
    input:
        bam = BAMDIR + "/{sample}.dedup.sorted.bam"
    output:
        bw = OUTDIR + "/1_bw/{sample}_raw_both.bw"
    log:
        OUTDIR + "/1_bw/{sample}.log"
    params:
        prefix = OUTDIR + "/1_bw/{sample}"
    shell:
        """
        ./scripts/bam2bw.sh {input.bam} {params.prefix} &> {log}
        """
