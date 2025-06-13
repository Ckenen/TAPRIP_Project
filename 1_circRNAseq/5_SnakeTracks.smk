#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
# SAMPLES = ["20210720_circRNAseq_oeHNRNPK_Rep1"]
BAMDIR1 = "results/04_CIRCexplorer/03_tophat2_filtered"
BAMDIR2 = "results/04_CIRCexplorer/05_tophat2_fusion"
OUTDIR = "results/05_tracks"

rule all:
    input:
        expand(OUTDIR + "/01_bw_rna_seq/{sample}_raw_both.bw", sample=SAMPLES),
        expand(OUTDIR + "/02_bw_fusion/{sample}_raw_both.bw", sample=SAMPLES),
        expand(OUTDIR + "/03_bw_fusion_clean/{sample}_raw_both.bw", sample=SAMPLES),


rule bam2bw_rna_seq:
    input:
        bam = BAMDIR1 + "/{sample}.bam"
    output:
        bw = OUTDIR + "/01_bw_rna_seq/{sample}_raw_both.bw"
    log:
        OUTDIR + "/01_bw_rna_seq/{sample}.log"
    params:
        prefix = OUTDIR + "/01_bw_rna_seq/{sample}"
    shell:
        """
        ./scripts/tracks/bam2bw_stranded_pe.sh {input.bam} {params.prefix} &> {log}
        """

rule bam2bw_fusion:
    input:
        bam = BAMDIR2 + "/{sample}/accepted_hits.bam"
    output:
        bw = OUTDIR + "/02_bw_fusion/{sample}_raw_both.bw"
    log:
        OUTDIR + "/02_bw_fusion/{sample}.log"
    params:
        prefix = OUTDIR + "/02_bw_fusion/{sample}"
    shell:
        """
        ./scripts/tracks/bam2bw_stranded_pe.sh {input.bam} {params.prefix} &> {log}
        """

rule bam2bw_fusion_clean:
    input:
        bam = BAMDIR2 + "/{sample}/accepted_hits.bam"
    output:
        bam1 = temp(OUTDIR + "/03_bw_fusion_clean/{sample}_clean.bam"),
        bam2 = OUTDIR + "/03_bw_fusion_clean/{sample}_clean.sorted.bam",
        bw = OUTDIR + "/03_bw_fusion_clean/{sample}_raw_both.bw"
    log:
        OUTDIR + "/03_bw_fusion_clean/{sample}.log"
    params:
        prefix = OUTDIR + "/03_bw_fusion_clean/{sample}"
    threads:
        4
    shell:
        """(
        ./scripts/mapping/filter_fusion_mapped_bam.py {input.bam} {output.bam1}
        samtools sort -@ {threads} -o {output.bam2} {output.bam1}
        samtools index -@ {threads} {output.bam2}
        ./scripts/tracks/bam2bw_stranded_pe.sh {output.bam2} {params.prefix} ) &> {log}
        """