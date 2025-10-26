#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
OUTDIR = "results/5_tracks"

rule all:
    input:
        expand(OUTDIR + "/1_bw_rna_seq/{sample}_raw_both.bw", sample=SAMPLES),
        expand(OUTDIR + "/2_bw_fusion/{sample}_raw_both.bw", sample=SAMPLES),
        expand(OUTDIR + "/3_bw_fusion_clean/{sample}_raw_both.bw", sample=SAMPLES),


rule bam2bw_rna_seq:
    input:
        bam = "results/2_CIRCexplorer/3_tophat2_filtered/{sample}.bam"
    output:
        bw = OUTDIR + "/1_bw_rna_seq/{sample}_raw_both.bw"
    log:
        OUTDIR + "/1_bw_rna_seq/{sample}.log"
    params:
        prefix = OUTDIR + "/1_bw_rna_seq/{sample}"
    shell:
        """
        ./scripts/bam2bw_stranded_pe.sh {input.bam} {params.prefix} &> {log}
        """

rule bam2bw_fusion:
    input:
        bam = "results/2_CIRCexplorer/5_tophat2_fusion/{sample}/accepted_hits.bam"
    output:
        bw = OUTDIR + "/2_bw_fusion/{sample}_raw_both.bw"
    log:
        OUTDIR + "/2_bw_fusion/{sample}.log"
    params:
        prefix = OUTDIR + "/2_bw_fusion/{sample}"
    shell:
        """
        ./scripts/bam2bw_stranded_pe.sh {input.bam} {params.prefix} &> {log}
        """

rule bam2bw_fusion_clean:
    input:
        bam = "results/2_CIRCexplorer/5_tophat2_fusion/{sample}/accepted_hits.bam"
    output:
        bam1 = temp(OUTDIR + "/3_bw_fusion_clean/{sample}_clean.bam"),
        bam2 = OUTDIR + "/3_bw_fusion_clean/{sample}_clean.sorted.bam",
        bw = OUTDIR + "/3_bw_fusion_clean/{sample}_raw_both.bw"
    log:
        OUTDIR + "/3_bw_fusion_clean/{sample}.log"
    params:
        prefix = OUTDIR + "/3_bw_fusion_clean/{sample}"
    threads:
        4
    shell:
        """(
        ./scripts/filter_fusion_mapped_bam.py {input.bam} {output.bam1}
        samtools sort -@ {threads} -o {output.bam2} {output.bam1}
        samtools index -@ {threads} {output.bam2}
        ./scripts/bam2bw_stranded_pe.sh {output.bam2} {params.prefix} ) &> {log}
        """