#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
# SAMPLES = ["20210720_circRNAseq_oeHNRNPK_Rep1"]
SAMPLES2 = list(filter(lambda s: "_circRNAseq_" in s, SAMPLES))
FQDIR = "results/01_prepare/03_bowtie2_mapping_ribo"
OUTDIR = "results/04_CIRCexplorer"

rule all:
    input:
        OUTDIR + "/01_genome_index/bowtie_index",
        OUTDIR + "/01_genome_index/bowtie2_index",
        OUTDIR + "/01_genome_index/tophat2_transcriptome_index",
        expand(OUTDIR + "/02_tophat2_mapped/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/03_tophat2_filtered/{sample}.bam", sample=SAMPLES), # optional
        expand(OUTDIR + "/04_tophat2_unmapped/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/05_tophat2_fusion/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/06_circ_parse/{sample}.bed", sample=SAMPLES),
        expand(OUTDIR + "/07_circ_annotate/{sample}.txt", sample=SAMPLES),
        expand(OUTDIR + "/08_circ_assemble/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/09_circ_denovo/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/10_circ_quantify/{sample}.txt", sample=SAMPLES),
        expand(OUTDIR + "/11_circ_quantify_vs_rnaseq/{sample}.txt", sample=SAMPLES2),


rule build_bowtie_index:
    input:
        fasta = config["GENOME_FASTA"]
    output:
        out = directory(OUTDIR + "/01_genome_index/bowtie_index")
    log:
        log = OUTDIR + "/01_genome_index/bowtie_index.log"
    threads:
        12
    conda:
        "tophat2"
    shell:
        """
        mkdir {output.out}
        bowtie-build --threads {threads} {input.fasta} {output}/ref &> {log}
        """

rule build_bowtie2_index:
    input:
        fasta = config["GENOME_FASTA"]
    output:
        out = directory(OUTDIR + "/01_genome_index/bowtie2_index")
    log:
        log = OUTDIR + "/01_genome_index/bowtie2_index.log"
    threads:
        12
    conda:
        "tophat2"
    shell:
        """
        mkdir {output.out}
        bowtie2-build --threads {threads} {input.fasta} {output}/ref &> {log}
        """

rule tophat2_transcriptome_index:
    input:
        gtf = config["GENOME_GTF"],
        bt2 = rules.build_bowtie2_index.output
    output:
        out = directory(OUTDIR + "/01_genome_index/tophat2_transcriptome_index")
    log:
        OUTDIR + "/01_genome_index/tophat2_transcriptome_index.log"
    threads:
        12
    conda:
        "tophat2"
    shell:
        """
        tophat -p {threads} -G {input.gtf} --transcriptome-index={output}/ref {input.bt2}/ref &> {log}
        """

rule tophat2: # The version of TopHat is v2.1.0, v2.1.1 will occur error in tophat-fusion running.
    input:
        fq1 = FQDIR + "/{sample}.unmapped.fastq.1.gz",
        fq2 = FQDIR + "/{sample}.unmapped.fastq.2.gz",
        bt2 = rules.build_bowtie2_index.output.out,
        th2 = rules.tophat2_transcriptome_index.output.out
    output:
        out = directory(OUTDIR + "/02_tophat2_mapped/{sample}")
    log:
        log = OUTDIR + "/02_tophat2_mapped/{sample}.log"
    threads:
        24
    conda:
        "tophat2"
    shell:
        """
        tophat2 -a 6 --microexon-search -m 2 -g 1 -p {threads} --library-type fr-firststrand \
            --no-discordant --no-mixed --transcriptome-index={input.th2}/ref \
            -o {output.out} {input.bt2}/ref {input.fq1} {input.fq2} &> {log}
        samtools index -@ {threads} {output}/accepted_hits.bam
        """

rule filter_bam:
    input:
        rules.tophat2.output.out
    output:
        bam = OUTDIR + "/03_tophat2_filtered/{sample}.bam"
    log:
        OUTDIR + "/03_tophat2_filtered/{sample}.log"
    threads:
        4
    conda:
        "tophat2"
    shell:
        """
        samtools view -@ {threads} -f 2 -F 2316 -q 30 -b -o {output.bam} {input}/accepted_hits.bam &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule get_unmapped_fq:
    input:
        rules.tophat2.output.out
    output:
        directory(OUTDIR + "/04_tophat2_unmapped/{sample}")
    threads:
        8
    conda:
        "tophat2"
    shell:
        """
        mkdir {output}
        samtools sort -@ {threads} -n {input}/unmapped.bam > {output}/sorted_by_name.bam
        bamToFastq -i {output}/sorted_by_name.bam -fq {output}/mate1.fastq -fq2 {output}/mate2.fastq
        pigz -p {threads} {output}/mate1.fastq {output}/mate2.fastq
        rm {output}/sorted_by_name.bam
        """

rule tophat2_fusion:
    input:
        fqs = rules.get_unmapped_fq.output,
        idx = rules.build_bowtie_index.output.out
    output:
        out = directory(OUTDIR + "/05_tophat2_fusion/{sample}")
    log:
        log = OUTDIR + "/05_tophat2_fusion/{sample}.log"
    threads:
        24
    conda:
        "tophat2"
    shell:
        """
        tophat2 -o {output.out} -p {threads} --fusion-search --keep-fasta-order \
            --bowtie1 --no-coverage-search {input.idx}/ref \
            {input.fqs}/mate1.fastq.gz {input.fqs}/mate2.fastq.gz &> {log}
        """

rule circ_parse:
    input:
        rules.tophat2_fusion.output.out
    output:
        bed = OUTDIR + "/06_circ_parse/{sample}.bed"
    log:
        OUTDIR + "/06_circ_parse/{sample}.log"
    conda:
        "CIRCexplorer2"
    shell:
        """
        CIRCexplorer2 parse -b {output} -t TopHat-Fusion {input}/accepted_hits.bam &> {log}
        """

rule circ_annotate:
    input:
        fa = config["GENOME_FASTA"],
        txt = config["GENEPRED"],
        bed = rules.circ_parse.output.bed
    output:
        txt = OUTDIR + "/07_circ_annotate/{sample}.txt"
    log:
        OUTDIR + "/07_circ_annotate/{sample}.log"
    conda:
        "CIRCexplorer2"
    shell:
        """
        CIRCexplorer2 annotate -r {input.txt} -g {input.fa} -b {input.bed} -o {output.txt} &> {log}
        """

rule circ_assemble:
    input:
        txt = config["GENEPRED"],
        bam = rules.tophat2_fusion.output.out
    output:
        directory(OUTDIR + "/08_circ_assemble/{sample}")
    log:
        OUTDIR + "/08_circ_assemble/{sample}.log"
    threads:
        12
    conda:
        "CIRCexplorer2"
    shell:
        """
        CIRCexplorer2 assemble -p {threads} -r {input.txt} -m {input.bam} -o {output} &> {log}
        """

# 如果想检测alternative splicing的话，需要poly(A)+组分作为input，现在只检测alternative back-splicing。

rule circ_denovo:
    input:
        genome = config["GENOME_FASTA"],
        ref = config["GENEPRED"],
        junc = rules.circ_parse.output.bed,
        cuff = rules.circ_assemble.output
    output:
        d1 = directory(OUTDIR + "/09_circ_denovo/{sample}"),
        d2 = directory(OUTDIR + "/09_circ_denovo/{sample}.abs")
    log:
        OUTDIR + "/09_circ_denovo/{sample}.log"
    conda:
        "CIRCexplorer2"
    shell:
        """
        CIRCexplorer2 denovo --abs={output.d2} -r {input.ref} -g {input.genome} -b {input.junc} -d {input.cuff} -o {output.d1} &> {log}
        """

# rule txt2bed:
#     input:
#         OUTDIR + "/circ_explorer2/denovo/{group}/{sample}"
#     output:
#         bed = OUTDIR + "/circRNA/{group}/{sample}.bed.gz"
#     shell:
#         """
#         awk 'BEGIN{{FS="\\t";OFS="\\t"}}{{print $1,$2,$3,$14"/"$15,$13,$6,$7,$8,$9,$10,$11,$12}}' \
#             {input}/circularRNA_full.txt | sort -k1,1 -k2,2n | bgzip -c > {output.bed}
#         tabix -p bed {output.bed}
#         """

rule quantify:
    input:
        circ = rules.circ_denovo.output.d1,
        bam = rules.tophat2.output.out,
        gene = config["GENEPRED"]
    output:
        txt = OUTDIR + "/10_circ_quantify/{sample}.txt"
    log:
        OUTDIR + "/10_circ_quantify/{sample}.log"
    conda:
        "CIRCexplorer2"
    shell:
        """
        circ_quant -c {input.circ}/circularRNA_full.txt -b {input.bam}/accepted_hits.bam -r {input.gene} -o {output.txt} &> {log}
        """

def get_rnaseq_bam(sample):
    sample = sample.replace("circRNAseq", "RNAseq")
    return OUTDIR + "/02_tophat2_mapped/%s" % sample

rule quantify_vs_rnaseq:
    input:
        circ = rules.circ_denovo.output.d1,
        bam = lambda wildcards: get_rnaseq_bam(wildcards.sample),
        gene = config["GENEPRED"]
    output:
        txt = OUTDIR + "/11_circ_quantify_vs_rnaseq/{sample}.txt"
    log:
        OUTDIR + "/11_circ_quantify_vs_rnaseq/{sample}.log"
    conda:
        "CIRCexplorer2"
    shell:
        """
        circ_quant -c {input.circ}/circularRNA_full.txt -b {input.bam}/accepted_hits.bam -r {input.gene} -o {output.txt} &> {log}
        """
