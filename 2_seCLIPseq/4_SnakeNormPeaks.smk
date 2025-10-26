#!/usr/bin/env runsnakemake
configfile: "config.yaml"
DATA = {
    "20241028_293T_HNRNPK_CLIP_Rep1": ["20241028_293T_HNRNPK_CLIP_IP_Rep1", "20241028_293T_HNRNPK_CLIP_Input_Rep1"],
    "20241028_293T_HNRNPK_CLIP_Rep2": ["20241028_293T_HNRNPK_CLIP_IP_Rep2", "20241028_293T_HNRNPK_CLIP_Input_Rep2"],
    "20241028_293T_RBM12_CLIP_Rep1": ["20241028_293T_RBM12_CLIP_IP_Rep1", "20241028_293T_RBM12_CLIP_Input_Rep1"],
    "20241028_293T_RBM12_CLIP_Rep2": ["20241028_293T_RBM12_CLIP_IP_Rep2", "20241028_293T_RBM12_CLIP_Input_Rep2"],
}
NAMES = list(sorted(DATA.keys()))
BAMDIR = "results/1_pipeline/9_dedup"
PEAKDIR = "results/2_peaks/1_clipper"
OUTDIR = "results/4_norm_peaks"

rule all:
    input:
        expand(OUTDIR + "/1_peaks/{name}", name=NAMES),
        expand(OUTDIR + "/2_sequences/{name}.fasta", name=NAMES),
        expand(OUTDIR + "/3_homer2/{name}", name=NAMES),

rule norm_peaks:
    input:
        bam_ip = lambda wildcards: "%s/%s.dedup.sorted.bam" % (BAMDIR, DATA[wildcards.name][0]),
        bam_input = lambda wildcards: "%s/%s.dedup.sorted.bam" % (BAMDIR, DATA[wildcards.name][1]),
        peak_ip = lambda wildcards: "%s/%s.bed.gz" % (PEAKDIR, DATA[wildcards.name][0]),
    output:
        directory(OUTDIR + "/1_peaks/{name}")
    log:
        OUTDIR + "/1_peaks/{name}.log"
    shell:
        """
        sh ./scripts/norm_peaks.sh {input} {output} &> {log}
        """

rule get_peak_fasta:
    input:
        beddir = rules.norm_peaks.output,
        fa = config["GENOME_FASTA"]
    output:
        bed = OUTDIR + "/2_sequences/{name}.bed",
        fa = OUTDIR + "/2_sequences/{name}.fasta"
    shell:
        """
        cat {input.beddir}/peakClusters.normed.compressed.bed | awk '$5>0.5' > {output.bed}
        bedtools getfasta -s -nameOnly -fi {input.fa} -bed {output.bed} > {output.fa}
        """

rule find_motifs:
    input:
        fa = rules.get_peak_fasta.output.fa
    output:
        out = directory(OUTDIR + "/3_homer2/{name}")
    log:
        OUTDIR + "/3_homer2/{name}.log"
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