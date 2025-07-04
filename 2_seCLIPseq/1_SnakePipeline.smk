#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
OUTDIR = "results/01_pipeline"

## seCLIP-seq pipeline.

rule all:
    input:
        # expand(OUTDIR + "/01_umitools/{sample}.fastq.gz", sample=SAMPLES),
        # expand(OUTDIR + "/02_cutadapt/{sample}.fastq.gz", sample=SAMPLES),
        # expand(OUTDIR + "/03_cutadapt2/{sample}.fastq.gz", sample=SAMPLES),
        # expand(OUTDIR + "/04_sorted_read_name/{sample}.fastq.gz", sample=SAMPLES),
        # expand(OUTDIR + "/05_bowtie2_mapped_ribo/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/06_star_mapped_repetitive/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/07_star_mapped_genome/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/08_sorted_bam/{sample}.byPos.bam", sample=SAMPLES),
        expand(OUTDIR + "/08_sorted_bam/{sample}.byPos.flagstat", sample=SAMPLES),
        expand(OUTDIR + "/09_dedup/{sample}.dedup.sorted.bam", sample=SAMPLES),
        expand(OUTDIR + "/09_dedup/{sample}.dedup.sorted.flagstat", sample=SAMPLES),

def get_fastq(sample):
    path = "data/datasets/%s_R1.fastq.gz" % sample
    if not os.path.exists(path):
        path = "data/public/GSE180686/%s.fastq.gz" % sample
    return path

rule umi_tools:
    input:
        fq = lambda wildcards: get_fastq(wildcards.sample)
    output:
        fq = temp(OUTDIR + "/01_umitools/{sample}.fastq.gz")
    log:
        OUTDIR + "/01_umitools/{sample}.log"
    shell:
        """
        umi_tools extract --random-seed 1 --bc-pattern NNNNNNNNNN \
            --stdin {input.fq} --stdout {output.fq} --log {log}
        """

rule cutadapt:
    input:
        fq = rules.umi_tools.output.fq
    output:
        fq = temp(OUTDIR + "/02_cutadapt/{sample}.fastq.gz")
    log:
        OUTDIR + "/02_cutadapt/{sample}.log"
    conda:
        "cutadapt"
    threads:
        8
    shell:
        """
        cutadapt -O 1 -j {threads} --match-read-wildcards \
            --times 1 -e 0.1 --quality-cutoff 6 -m 18 -o {output.fq} \
            -a AGATCGGAAGAGCAC \
            -a GATCGGAAGAGCACA \
            -a ATCGGAAGAGCACAC \
            -a TCGGAAGAGCACACG \
            -a CGGAAGAGCACACGT \
            -a GGAAGAGCACACGTC \
            -a GAAGAGCACACGTCT \
            -a AAGAGCACACGTCTG \
            -a AGAGCACACGTCTGA \
            -a GAGCACACGTCTGAA \
            -a AGCACACGTCTGAAC \
            -a GCACACGTCTGAACT \
            -a CACACGTCTGAACTC \
            -a ACACGTCTGAACTCC \
            -a CACGTCTGAACTCCA \
            -a ACGTCTGAACTCCAG \
            -a CGTCTGAACTCCAGT \
            -a GTCTGAACTCCAGTC \
            -a TCTGAACTCCAGTCA \
            -a CTGAACTCCAGTCAC \
            {input.fq} &> {log}
        """

rule cutadapt2:
    input:
        fq = rules.cutadapt.output.fq
    output:
        fq = temp(OUTDIR + "/03_cutadapt2/{sample}.fastq.gz")
    log:
        OUTDIR + "/03_cutadapt2/{sample}.log"
    conda:
        "cutadapt"
    threads:
        8
    shell:
        """
        cutadapt -O 5 -j {threads} --match-read-wildcards --times 1 \
            -e 0.1 --quality-cutoff 6 -m 18 -o {output.fq} \
            -a AGATCGGAAGAGCAC \
            -a GATCGGAAGAGCACA \
            -a ATCGGAAGAGCACAC \
            -a TCGGAAGAGCACACG \
            -a CGGAAGAGCACACGT \
            -a GGAAGAGCACACGTC \
            -a GAAGAGCACACGTCT \
            -a AAGAGCACACGTCTG \
            -a AGAGCACACGTCTGA \
            -a GAGCACACGTCTGAA \
            -a AGCACACGTCTGAAC \
            -a GCACACGTCTGAACT \
            -a CACACGTCTGAACTC \
            -a ACACGTCTGAACTCC \
            -a CACGTCTGAACTCCA \
            -a ACGTCTGAACTCCAG \
            -a CGTCTGAACTCCAGT \
            -a GTCTGAACTCCAGTC \
            -a TCTGAACTCCAGTCA \
            -a CTGAACTCCAGTCAC \
            {input.fq} &> {log}
        """

rule sort_fastq:
    input:
        fq = rules.cutadapt2.output.fq
    output:
        fq = temp(OUTDIR + "/04_sorted_read_name/{sample}.fastq.gz")
    shell:
        """
        sort_fastq_name.sh {input.fq} | gzip -c > {output.fq}
        """

rule bowtie2_mapping_ribo:
    input:
        fq = rules.sort_fastq.output.fq,
        idx = config["BOWTIE2_INDEX_RIBO"]
    output:
        bam = OUTDIR + "/05_bowtie2_mapped_ribo/{sample}.bam",
        fq = OUTDIR + "/05_bowtie2_mapped_ribo/{sample}.fastq.gz",
    log:
        OUTDIR + "/05_bowtie2_mapped_ribo/{sample}.log"
    params:
        prefix = OUTDIR + "/05_bowtie2_mapped_ribo/{sample}"
    threads:
        12
    conda:
        "bowtie2"
    shell:
        """(
        bowtie2 -p {threads} --no-unal -t --un-gz {output.fq} -x {input.idx}/ref -U {input.fq} \
            | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -o {output.bam} - ) &> {log}
        """

rule star_mapping_repetitive:
    input:
        fq = rules.bowtie2_mapping_ribo.output.fq,
        idx = config["STAR_INDEX_REP"]
    output:
        out = directory(OUTDIR + "/06_star_mapped_repetitive/{sample}")
    log:
        OUTDIR + "/06_star_mapped_repetitive/{sample}.log"
    threads:
        12
    conda:
        "star"
    shell:
        """
        mkdir -p {output}
        STAR --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {input.idx} \
            --alignEndsType EndToEnd \
            --outSAMunmapped Within \
            --outFilterMultimapNmax 30 \
            --outFilterMultimapScoreRange 1 \
            --outFileNamePrefix {output}/{wildcards.sample}. \
            --outSAMattributes All \
            --outSAMtype BAM Unsorted \
            --outFilterType BySJout \
            --outReadsUnmapped Fastx \
            --outFilterScoreMin 10 \
            --outSAMattrRGline ID:foo \
            --outStd Log \
            --outBAMcompression 10 \
            --outSAMmode Full \
            --readFilesCommand zcat \
            --readFilesIn {input.fq} &> {log}
        """

rule star_mapping_genome:
    input:
        fqs = rules.star_mapping_repetitive.output.out,
        idx = config["STAR_INDEX_GENOME"]
    output:
        out = directory(OUTDIR + "/07_star_mapped_genome/{sample}")
    log:
        OUTDIR + "/07_star_mapped_genome/{sample}.log"
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
            --alignEndsType EndToEnd \
            --outSAMunmapped Within \
            --outFilterMultimapNmax 1 \
            --outFilterMultimapScoreRange 1 \
            --outFileNamePrefix {output}/{wildcards.sample}. \
            --outSAMattributes All \
            --outSAMtype BAM Unsorted \
            --outFilterType BySJout \
            --outReadsUnmapped Fastx \
            --outFilterScoreMin 10 \
            --outSAMattrRGline ID:foo \
            --outStd Log \
            --outBAMcompression 10 \
            --outSAMmode Full \
            --readFilesIn \
            {input.fqs}/{wildcards.sample}.Unmapped.out.mate1 &> {log}
        """

rule sort_bam:
    input:
        bamdir = rules.star_mapping_genome.output.out
    output:
        bam1 = temp(OUTDIR + "/08_sorted_bam/{sample}.byName.bam"),
        bam2 = OUTDIR + "/08_sorted_bam/{sample}.byPos.bam"
    threads:
        4
    shell:
        """
        samtools sort -@ {threads} -n -o {output.bam1} {input.bamdir}/{wildcards.sample}.Aligned.out.bam
        samtools sort -@ {threads} -o {output.bam2} {output.bam1}
        samtools index -@ {threads} {output.bam2}
        """

rule dedup:
    input:
        bam = rules.sort_bam.output.bam2
    output:
        bam1 = temp(OUTDIR + "/09_dedup/{sample}.dedup.bam"),
        bam2 = OUTDIR + "/09_dedup/{sample}.dedup.sorted.bam"
    log:
        OUTDIR + "/09_dedup/{sample}.dedup.log"
    params:
        prefix = OUTDIR + "/09_dedup/{sample}.dedup.stats"
    threads:
        4
    shell:
        """(
        umi_tools dedup --random-seed 1 -I {input.bam} --method unique \
            --output-stats {params.prefix} -S {output.bam1}
        samtools sort -@ {threads} -o {output.bam2} {output.bam1} 
        samtools index -@ {threads} {output.bam2} ) &> {log}
        """

rule flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    threads:
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """