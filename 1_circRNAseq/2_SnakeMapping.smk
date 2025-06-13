#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
# SAMPLES = ["20210720_RNAseq_oeHNRNPK_Rep2"]
FQDIR = "results/01_prepare/03_bowtie2_mapping_ribo"
OUTDIR = "results/02_mapping"

rule all:
    input:
        OUTDIR + "/01_star_index_genome/genome_star_index",
        expand(OUTDIR + "/02_star_mapping/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/03_filtered/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/03_filtered/{sample}.flagstat", sample=SAMPLES),
        expand(OUTDIR + "/03_filtered/{sample}.infer_experiment.txt", sample=SAMPLES),
        expand(OUTDIR + "/04_rmdup/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/04_rmdup/{sample}.flagstat", sample=SAMPLES),


rule star_build_genome:
    input:
        fa = config["GENOME_FASTA"],
        gtf = config["GENOME_GTF"]
    output:
        idx = directory(OUTDIR + "/01_star_index_genome/genome_star_index")
    log:
        OUTDIR + "/01_star_index_genome/genome_star_index.log"
    threads:
        12
    conda:
        "star"
    shell:
        """
        mkdir -p {output.idx}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output.idx} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} &> {log}
        """

# STAR --runMode alignReads --genomeDir results/mapping/star/index --genomeLoad Remove

rule star_mapping:
    input:
        fq1 = FQDIR + "/{sample}.unmapped.fastq.1.gz",
        fq2 = FQDIR + "/{sample}.unmapped.fastq.2.gz",
        idx = rules.star_build_genome.output.idx
    output:
        out = directory(OUTDIR + "/02_star_mapping/{sample}")
    log:
        OUTDIR + "/02_star_mapping/{sample}.log"
    threads:
        12
    conda:
        "star"
    shell:
        """
        mkdir {output}
        STAR --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {input.idx} \
            --genomeLoad LoadAndKeep \
            --outFileNamePrefix {output}/{wildcards.sample}. \
            --outReadsUnmapped Fastx \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes All \
            --limitBAMsortRAM 64000000000 \
            --readFilesCommand zcat \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        samtools index -@ {threads} {output}/{wildcards.sample}.Aligned.sortedByCoord.out.bam
        pigz -p {threads} {output}/*.mate*
        """

rule filter_bam:
    input:
        rules.star_mapping.output.out
    output:
        bam = OUTDIR + "/03_filtered/{sample}.bam"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} --expr 'rname =~ "^chr([0-9]+|[XY])$"' \
            -d 'NH:1' -q 30 -f 2 -F 2308 \
            -o {output.bam} {input}/{wildcards.sample}.Aligned.sortedByCoord.out.bam
        samtools index -@ {threads} {output.bam}
        """

rule infer_experiment:
    input:
        bam = rules.filter_bam.output.bam,
        bed = config["TRANSCRIPT_BED"]
    output:
        txt = OUTDIR + "/03_filtered/{sample}.infer_experiment.txt"
    shell:
        """
        infer_experiment.py -i {input.bam} -r {input.bed} > {output.txt} 2> /dev/null
        """

rule rmdup: 
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = OUTDIR + "/04_rmdup/{sample}.bam"
    log:
        OUTDIR + "/04_rmdup/{sample}.log"
    threads:
        4
    shell:
        """
        sambamba markdup -t {threads} -r {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    shell:
        """
        samtools flagstat {input.bam} > {output.txt}
        """
