#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
BAMDIR = "results/1_pipeline/9_dedup"
OUTDIR = "results/2_peaks"


rule all:
    input:
        expand(OUTDIR + "/1_clipper/{sample}.bed.gz", sample=SAMPLES),

# docker

rule clipper:
    input:
        bam = BAMDIR + "/{sample}.dedup.sorted.bam",
        bai = BAMDIR + "/{sample}.dedup.sorted.bam.bai",
    output:
        tmp = temp(directory(OUTDIR + "/1_clipper/{sample}.d")),
        bed = OUTDIR + "/1_clipper/{sample}.bed.gz"
    log:
        OUTDIR + "/1_clipper/{sample}.log"
    threads:
        24
    shell:
        """(
        mkdir -p {output.tmp}
        cp {input.bam} {output.tmp}/{wildcards.sample}.bam
        cp {input.bai} {output.tmp}/{wildcards.sample}.bam.bai
        docker run --rm -v `realpath {output.tmp}`:/data brianyee/clipper \
            clipper --bam /data/{wildcards.sample}.bam \
            --species GRCh38_v29e \
            --outfile /data/{wildcards.sample}.bed \
            --processors={threads} 
        sort -k1,1 -k2,2n {output.tmp}/{wildcards.sample}.bed | bgzip -c > {output.bed}
        tabix -p bed {output.bed} ) &> {log}
        """

# conda 

# rule clipper:
#     input:
#         bam = BAMDIR + "/{sample}.dedup.sorted.bam",
#         bai = BAMDIR + "/{sample}.dedup.sorted.bam.bai",
#     output:
#         tmp = temp(directory(OUTDIR + "/1_clipper/{sample}.d")),
#         bed = OUTDIR + "/1_clipper/{sample}.bed.gz"
#     log:
#         OUTDIR + "/1_clipper/{sample}.log"
#     threads:
#         24
#     conda:
#         "clipper3"
#     shell:
#         """(
#         clipper --bam {input.bam} --species GRCh38_v29e --outfile {output.bed1} --processors={threads} 
#         sort -k1,1 -k2,2n {output.bed1} | bgzip -c > {output.bed2}
#         tabix -p bed {output.bed2} ) &> {log}
#         """
