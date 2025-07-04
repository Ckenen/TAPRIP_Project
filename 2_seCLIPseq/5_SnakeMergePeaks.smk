#!/usr/bin/env runsnakemake
configfile: "config.yaml"
DATA = {
    "20241028_293T_HNRNPK_CLIP": [
        "20241028_293T_HNRNPK_CLIP_IP_Rep1",
        "20241028_293T_HNRNPK_CLIP_IP_Rep2",
        "20241028_293T_HNRNPK_CLIP_Input_Rep1",
        "20241028_293T_HNRNPK_CLIP_Input_Rep2",
    ],
    "20241028_293T_RBM12_CLIP": [
        "20241028_293T_RBM12_CLIP_IP_Rep1", 
        "20241028_293T_RBM12_CLIP_IP_Rep2", 
        "20241028_293T_RBM12_CLIP_Input_Rep1", 
        "20241028_293T_RBM12_CLIP_Input_Rep2",
    ],
    "20241125_seCLIPseq_293T_Ctrl": [
        "20241125_seCLIPseq_293T_Ctrl_IP_Rep1", 
        "20241125_seCLIPseq_293T_Ctrl_IP_Rep2", 
        "20241125_seCLIPseq_293T_Ctrl_Input_Rep1", 
        "20241125_seCLIPseq_293T_Ctrl_Input_Rep2",
    ],
    "20241125_seCLIPseq_293T_XIST_A": [
        "20241125_seCLIPseq_293T_XIST_A_IP_Rep1", 
        "20241125_seCLIPseq_293T_XIST_A_IP_Rep2", 
        "20241125_seCLIPseq_293T_XIST_A_Input_Rep1", 
        "20241125_seCLIPseq_293T_XIST_A_Input_Rep2",
    ],
    "20241125_seCLIPseq_293T_XIST_E": [
        "20241125_seCLIPseq_293T_XIST_E_IP_Rep1", 
        "20241125_seCLIPseq_293T_XIST_E_IP_Rep2", 
        "20241125_seCLIPseq_293T_XIST_E_Input_Rep1", 
        "20241125_seCLIPseq_293T_XIST_E_Input_Rep2",
    ],
    "20241128_seCLIPseq_293T_RBM12": [
        "20241128_seCLIPseq_293T_RBM12_IP_Rep1",
        "20241128_seCLIPseq_293T_RBM12_IP_Rep2",
        "20241128_seCLIPseq_293T_RBM12_Input_Rep1",
        "20241128_seCLIPseq_293T_RBM12_Input_Rep2",
    ],
}
NAMES = list(sorted(DATA.keys()))
BAMDIR = "results/01_pipeline/09_dedup"
PEAKDIR = "results/02_peaks/01_clipper"
OUTDIR = "results/05_merged_rep_peaks"

rule all:
    input:
        expand(OUTDIR + "/01_peaks/{name}", name=NAMES),
        expand(OUTDIR + "/02_sequences/{name}.fasta", name=NAMES),
        expand(OUTDIR + "/03_homer2/{name}", name=NAMES),

rule merge_rep_peaks:
    input:
        bam_ip_1 = lambda wildcards: "%s/%s.dedup.sorted.bam" % (BAMDIR, DATA[wildcards.name][0]),
        bam_ip_2 = lambda wildcards: "%s/%s.dedup.sorted.bam" % (BAMDIR, DATA[wildcards.name][1]),
        bam_input_1 = lambda wildcards: "%s/%s.dedup.sorted.bam" % (BAMDIR, DATA[wildcards.name][2]),
        bam_input_2 = lambda wildcards: "%s/%s.dedup.sorted.bam" % (BAMDIR, DATA[wildcards.name][3]),
        peak_1 = lambda wildcards: "%s/%s.bed.gz" % (PEAKDIR, DATA[wildcards.name][0]),
        peak_2 = lambda wildcards: "%s/%s.bed.gz" % (PEAKDIR, DATA[wildcards.name][1]),
    output:
        directory(OUTDIR + "/01_peaks/{name}")
    log:
        OUTDIR + "/01_peaks/{name}.log"
    shell:
        """
        sh ./scripts/merge_rep_peaks.sh {input} {output} &> {log}
        """

rule get_peak_fasta:
    input:
        bed = rules.merge_rep_peaks.output,
        fa = config["GENOME_FASTA"]
    output:
        fa = OUTDIR + "/02_sequences/{name}.fasta"
    shell:
        """
        bedtools getfasta -s -nameOnly -fi {input.fa} -bed {input.bed}/01v02.idr.out.bed > {output.fa}
        """

rule find_motifs:
    input:
        fa = rules.get_peak_fasta.output.fa
    output:
        out = directory(OUTDIR + "/03_homer2/{name}")
    log:
        OUTDIR + "/03_homer2/{name}.log"
    threads:
        8
    conda:
        "homer2"
    shell:
        """
        findMotifs.pl {input.fa} fasta {output}/ -p {threads} -rna -norevopp -mask -len 6,7,8 &> {log}
        """