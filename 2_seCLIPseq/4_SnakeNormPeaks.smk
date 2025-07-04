#!/usr/bin/env runsnakemake
configfile: "config.yaml"
DATA = {
    "20241001_293T_RBM12_CLIP_Rep1": ["20241001_293T_RBM12_CLIP_IP_Rep1", "20241001_293T_RBM12_CLIP_Input_Rep1"],
    "20241001_293T_RBM12_CLIP_Rep2": ["20241001_293T_RBM12_CLIP_IP_Rep2", "20241001_293T_RBM12_CLIP_Input_Rep2"],
    "20241028_293T_HNRNPK_CLIP_Rep1": ["20241028_293T_HNRNPK_CLIP_IP_Rep1", "20241028_293T_HNRNPK_CLIP_Input_Rep1"],
    "20241028_293T_HNRNPK_CLIP_Rep2": ["20241028_293T_HNRNPK_CLIP_IP_Rep2", "20241028_293T_HNRNPK_CLIP_Input_Rep2"],
    "20241028_293T_RBM12_CLIP_Rep1": ["20241028_293T_RBM12_CLIP_IP_Rep1", "20241028_293T_RBM12_CLIP_Input_Rep1"],
    "20241028_293T_RBM12_CLIP_Rep2": ["20241028_293T_RBM12_CLIP_IP_Rep2", "20241028_293T_RBM12_CLIP_Input_Rep2"],
    "20241125_seCLIPseq_293T_Ctrl_Rep1": ["20241125_seCLIPseq_293T_Ctrl_IP_Rep1", "20241125_seCLIPseq_293T_Ctrl_Input_Rep1"],
    "20241125_seCLIPseq_293T_Ctrl_Rep2": ["20241125_seCLIPseq_293T_Ctrl_IP_Rep2", "20241125_seCLIPseq_293T_Ctrl_Input_Rep2"],
    "20241125_seCLIPseq_293T_XIST_A_Rep1": ["20241125_seCLIPseq_293T_XIST_A_IP_Rep1", "20241125_seCLIPseq_293T_XIST_A_Input_Rep1"],
    "20241125_seCLIPseq_293T_XIST_A_Rep2": ["20241125_seCLIPseq_293T_XIST_A_IP_Rep2", "20241125_seCLIPseq_293T_XIST_A_Input_Rep2"],
    "20241125_seCLIPseq_293T_XIST_E_Rep1": ["20241125_seCLIPseq_293T_XIST_E_IP_Rep1", "20241125_seCLIPseq_293T_XIST_E_Input_Rep1"],
    "20241125_seCLIPseq_293T_XIST_E_Rep2": ["20241125_seCLIPseq_293T_XIST_E_IP_Rep2", "20241125_seCLIPseq_293T_XIST_E_Input_Rep2"],
    "20241128_seCLIPseq_293T_RBM12_Rep1": ["20241128_seCLIPseq_293T_RBM12_IP_Rep1", "20241128_seCLIPseq_293T_RBM12_Input_Rep1"],
    "20241128_seCLIPseq_293T_RBM12_Rep2": ["20241128_seCLIPseq_293T_RBM12_IP_Rep2", "20241128_seCLIPseq_293T_RBM12_Input_Rep2"],
}
NAMES = list(sorted(DATA.keys()))
BAMDIR = "results/01_pipeline/09_dedup"
PEAKDIR = "results/02_peaks/01_clipper"
OUTDIR = "results/04_norm_peaks"

rule all:
    input:
        expand(OUTDIR + "/01_peaks/{name}", name=NAMES),
        expand(OUTDIR + "/02_sequences/{name}.fasta", name=NAMES),
        expand(OUTDIR + "/03_homer2/{name}", name=NAMES),

rule norm_peaks:
    input:
        bam_ip = lambda wildcards: "%s/%s.dedup.sorted.bam" % (BAMDIR, DATA[wildcards.name][0]),
        bam_input = lambda wildcards: "%s/%s.dedup.sorted.bam" % (BAMDIR, DATA[wildcards.name][1]),
        peak_ip = lambda wildcards: "%s/%s.bed.gz" % (PEAKDIR, DATA[wildcards.name][0]),
    output:
        directory(OUTDIR + "/01_peaks/{name}")
    log:
        OUTDIR + "/01_peaks/{name}.log"
    shell:
        """
        sh ./scripts/norm_peaks.sh {input} {output} &> {log}
        """

rule get_peak_fasta:
    input:
        beddir = rules.norm_peaks.output,
        fa = config["GENOME_FASTA"]
    output:
        bed = OUTDIR + "/02_sequences/{name}.bed",
        fa = OUTDIR + "/02_sequences/{name}.fasta"
    shell:
        """
        cat {input.beddir}/peakClusters.normed.compressed.bed | awk '$5>0.5' > {output.bed}
        bedtools getfasta -s -nameOnly -fi {input.fa} -bed {output.bed} > {output.fa}
        """

rule find_motifs:
    input:
        fa = rules.get_peak_fasta.output.fa
    output:
        out = directory(OUTDIR + "/03_homer2/{name}")
    log:
        OUTDIR + "/03_homer2/{name}.log"
    params:
        lens = lambda wildcards: "24,26,28,30,32,34,36,38,40" if "RBM12" in wildcards.name else "6,7,8"
    threads:
        8
    conda:
        "homer2"
    shell:
        """
        findMotifs.pl {input.fa} fasta {output}/ -p {threads} -rna -norevopp -mask -len {params.lens} &> {log}
        """