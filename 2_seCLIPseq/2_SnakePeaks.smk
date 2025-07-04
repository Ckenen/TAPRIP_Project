#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
# SAMPLES = list(filter(lambda s: "Input" not in s, SAMPLES))
BAMDIR = "results/01_pipeline/09_dedup"
OUTDIR = "results/02_peaks"


rule all:
    input:
        expand(OUTDIR + "/01_clipper/{sample}.bed.gz", sample=SAMPLES),

rule clipper:
    input:
        bam = BAMDIR + "/{sample}.dedup.sorted.bam"
    output:
        bed1 = temp(OUTDIR + "/01_clipper/{sample}.bed"),
        bed2 = OUTDIR + "/01_clipper/{sample}.bed.gz"
    log:
        OUTDIR + "/01_clipper/{sample}.log"
    conda:
        "clipper3"
    threads:
        24
    shell:
        """(
        clipper --bam {input.bam} --species GRCh38_v29e --outfile {output.bed1} --processors={threads} 
        sort -k1,1 -k2,2n {output.bed1} | bgzip -c > {output.bed2}
        tabix -p bed {output.bed2} ) &> {log}
        """
