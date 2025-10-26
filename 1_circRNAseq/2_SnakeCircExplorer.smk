#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
SAMPLES2 = list(filter(lambda s: "_circRNAseq_" in s, SAMPLES))
FQDIR = "results/1_prepare/3_bowtie2"
OUTDIR = "results/2_CIRCexplorer"


rule all:
    input:
        OUTDIR + "/1_genome_index/bowtie_index",
        OUTDIR + "/1_genome_index/bowtie2_index",
        OUTDIR + "/1_genome_index/tophat2_transcriptome_index",
        expand(OUTDIR + "/2_tophat2_mapped/{sample}", sample=SAMPLES),
        # expand(OUTDIR + "/3_tophat2_filtered/{sample}.bam", sample=SAMPLES), # optional
        expand(OUTDIR + "/4_tophat2_unmapped/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/5_tophat2_fusion/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/6_circ_parse/{sample}.bed", sample=SAMPLES),
        expand(OUTDIR + "/7_circ_annotate/{sample}.txt", sample=SAMPLES),
        expand(OUTDIR + "/8_circ_assemble/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/9_circ_denovo/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/10_circ_quantify/{sample}.txt", sample=SAMPLES),
        expand(OUTDIR + "/11_circ_quantify_vs_rnaseq/{sample}.txt", sample=SAMPLES2),

rule build_bowtie_index:
    input:
        fasta = config["FASTA"]
    output:
        out = directory(OUTDIR + "/1_genome_index/bowtie_index")
    log:
        log = OUTDIR + "/1_genome_index/bowtie_index.log"
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
        fasta = config["FASTA"]
    output:
        out = directory(OUTDIR + "/1_genome_index/bowtie2_index")
    log:
        log = OUTDIR + "/1_genome_index/bowtie2_index.log"
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
        gtf = config["GTF"],
        bt2 = rules.build_bowtie2_index.output
    output:
        out = directory(OUTDIR + "/1_genome_index/tophat2_transcriptome_index")
    log:
        OUTDIR + "/1_genome_index/tophat2_transcriptome_index.log"
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
        out = directory(OUTDIR + "/2_tophat2_mapped/{sample}")
    log:
        log = OUTDIR + "/2_tophat2_mapped/{sample}.log"
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
        bam = OUTDIR + "/3_tophat2_filtered/{sample}.bam"
    log:
        OUTDIR + "/3_tophat2_filtered/{sample}.log"
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
        directory(OUTDIR + "/4_tophat2_unmapped/{sample}")
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
        out = directory(OUTDIR + "/5_tophat2_fusion/{sample}")
    log:
        log = OUTDIR + "/5_tophat2_fusion/{sample}.log"
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
        bed = OUTDIR + "/6_circ_parse/{sample}.bed"
    log:
        OUTDIR + "/6_circ_parse/{sample}.log"
    conda:
        "CIRCexplorer2"
    shell:
        """
        CIRCexplorer2 parse -b {output} -t TopHat-Fusion {input}/accepted_hits.bam &> {log}
        """

rule circ_annotate:
    input:
        fa = config["FASTA"],
        txt = config["GENEPRED"],
        bed = rules.circ_parse.output.bed
    output:
        txt = OUTDIR + "/7_circ_annotate/{sample}.txt"
    log:
        OUTDIR + "/7_circ_annotate/{sample}.log"
    conda:
        "CIRCexplorer2"
    shell:
        """
        CIRCexplorer2 annotate -r {input.txt} \
            -g {input.fa} -b {input.bed} -o {output.txt} &> {log}
        """

rule circ_assemble:
    input:
        txt = config["GENEPRED"],
        bam = rules.tophat2_fusion.output.out
    output:
        directory(OUTDIR + "/8_circ_assemble/{sample}")
    log:
        OUTDIR + "/8_circ_assemble/{sample}.log"
    threads:
        12
    conda:
        "CIRCexplorer2"
    shell:
        """
        CIRCexplorer2 assemble -p {threads} \
            -r {input.txt} -m {input.bam} -o {output} &> {log}
        """

rule circ_denovo:
    input:
        genome = config["FASTA"],
        ref = config["GENEPRED"],
        junc = rules.circ_parse.output.bed,
        cuff = rules.circ_assemble.output
    output:
        d1 = directory(OUTDIR + "/9_circ_denovo/{sample}"),
        d2 = directory(OUTDIR + "/9_circ_denovo/{sample}.abs")
    log:
        OUTDIR + "/9_circ_denovo/{sample}.log"
    conda:
        "CIRCexplorer2"
    shell:
        """
        CIRCexplorer2 denovo --abs={output.d2} -r {input.ref} \
            -g {input.genome} -b {input.junc} \
            -d {input.cuff} -o {output.d1} &> {log}
        """

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
        circ_quant -c {input.circ}/circularRNA_full.txt \
            -b {input.bam}/accepted_hits.bam \
            -r {input.gene} -o {output.txt} &> {log}
        """

def get_rnaseq_bam(sample):
    sample = sample.replace("circRNAseq", "RNAseq")
    return OUTDIR + "/2_tophat2_mapped/%s" % sample

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
        circ_quant -c {input.circ}/circularRNA_full.txt \
            -b {input.bam}/accepted_hits.bam \
            -r {input.gene} -o {output.txt} &> {log}
        """
