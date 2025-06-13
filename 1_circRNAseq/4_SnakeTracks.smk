#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
BAMDIR1 = "results/02_mapping/04_rmdup"
OUTDIR = "results/04_tracks"

rule all:
    input:
        expand(OUTDIR + "/01_bw/{sample}_raw_both.bw", sample=SAMPLES),

rule bam2bw_rna_seq:
    input:
        bam = BAMDIR1 + "/{sample}.bam"
    output:
        bw = OUTDIR + "/01_bw/{sample}_raw_both.bw"
    log:
        OUTDIR + "/01_bw/{sample}.log"
    params:
        prefix = OUTDIR + "/01_bw/{sample}"
    shell:
        """
        ./scripts/tracks/bam2bw_stranded_pe.sh {input.bam} {params.prefix} &> {log}
        """