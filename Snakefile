import os
import yaml

configfile: 'CompProjectconfig.yaml' 

# Load metadata
single_samples = config["Single End"]
paired_samples = config["Paired End"]
samples = config["SRAs"]
reference_genome = config["Reference Genome"]

rule all:
    input:   
        expand("bigwig_files/pe/{sample}.bw", sample=paired_samples),
        expand("bigwig_files/se/{sample}.bw", sample=single_samples)

# --- FASTQ PREPARATION ---

rule fasterq_dump:    
    input:
        "initial_data/{sample}/{sample}.sra" 
    output:
        "data/fastq/pe/{sample}_1.fastq", 
        "data/fastq/pe/{sample}_2.fastq" 
    conda: "chipseq_env.yaml"
    shell:
        """
        mkdir -p data/fastq/pe
        fasterq-dump {input} -O data/fastq/pe --split-files
        """
        
rule fasterq_dump_single:
    input:
        "initial_data/{sample}/{sample}.sra" 
    output:
        "data/fastq/se/{sample}.fastq"
    conda: "chipseq_env.yaml" 
    shell:
        """
        mkdir -p data/fastq/se
        fasterq-dump {input} -O data/fastq/se
        """

# --- REFERENCE PREPARATION ---

rule download_reference_genome:
    output:
        genome="ref/P_falciparum3D7.fa",
        gff="ref/P_falciparum3D7_annotations.gff3"
    conda: "chipseq_env.yaml"
    shell:
        """
        mkdir -p ref 
        datasets download genome accession {reference_genome} --include gff3,genome --filename ref/ncbi_dataset.zip
        unzip -o ref/ncbi_dataset.zip -d ref/ncbi_dataset
        cp ref/ncbi_dataset/ncbi_dataset/data/*/*genomic.fna {output.genome}
        cp ref/ncbi_dataset/ncbi_dataset/data/*/*genomic.gff {output.gff}
        rm -rf ref/ncbi_dataset ref/ncbi_dataset.zip
        """

rule index_ref:
    input:
        genome="ref/P_falciparum3D7.fa"
    output:
        touch("ref/index.done")
    conda: "chipseq_env.yaml" 
    shell:
        """
        bwa index {input.genome}
        """

# --- TRIMMING ---

rule trimmomatic_se:
    input:
        "data/fastq/se/{sample}.fastq"
    output:
        "trimmed/se/{sample}.fastq"
    conda: "chipseq_env.yaml" 
    shell:
        """ 
        mkdir -p trimmed/se 
        trimmomatic SE -phred33 {input} {output} SLIDINGWINDOW:4:30 MINLEN:35
        """

rule trimmomatic_pe:
    input:
        r1="data/fastq/pe/{sample}_1.fastq",
        r2="data/fastq/pe/{sample}_2.fastq"
    output:
        r1_paired="trimmed/pe/{sample}_1_paired.fastq",
        r1_unpaired="trimmed/pe/{sample}_1_unpaired.fastq",
        r2_paired="trimmed/pe/{sample}_2_paired.fastq",
        r2_unpaired="trimmed/pe/{sample}_2_unpaired.fastq"
    conda: "chipseq_env.yaml"
    shell:
        """
        mkdir -p trimmed/pe
        trimmomatic PE -phred33 {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} SLIDINGWINDOW:4:30 MINLEN:35
        """

# --- MAPPING & FILTERING ---

rule bwa_mapping_se: 
    input:
        genome="ref/P_falciparum3D7.fa",
        fastq="trimmed/se/{sample}.fastq",
        donecheck="ref/index.done"
    output:
        "mapped_reads/se/{sample}.bam"
    conda: "chipseq_env.yaml"
    shell:
        """
        mkdir -p mapped_reads/se 
        bwa mem -M -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:lib1\\tPL:ILLUMINA' {input.genome} {input.fastq} | samtools view -bS - > {output}
        """

rule bwa_mapping_pe:
    input:
        genome="ref/P_falciparum3D7.fa",
        r1="trimmed/pe/{sample}_1_paired.fastq",
        r2="trimmed/pe/{sample}_2_paired.fastq",
        idx="ref/index.done"
    output:
        "mapped_reads/pe/{sample}.bam"
    conda: "chipseq_env.yaml"
    shell:
        """
        mkdir -p mapped_reads/pe
        bwa mem -M -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:lib1\\tPL:ILLUMINA' {input.genome} {input.r1} {input.r2} | samtools view -bS - > {output}
        """

rule filter_by_phred_score_se:
    input: "mapped_reads/se/{sample}.bam"
    output: "mapped_reads/se/{sample}.phred30.bam"
    conda: "chipseq_env.yaml"
    shell: "samtools view -b -q 30 {input} > {output}"

rule filter_by_phred_score_pe:
    input: "mapped_reads/pe/{sample}.bam"
    output: "mapped_reads/pe/{sample}.phred30.bam"
    conda: "chipseq_env.yaml"
    shell: "samtools view -b -q 30 {input} > {output}"

rule sort_bam_se:
    input: "mapped_reads/se/{sample}.phred30.bam"
    output: "mapped_reads/se/{sample}.sorted.bam"
    conda: "chipseq_env.yaml"
    shell: "samtools sort {input} -o {output}"

rule sort_bam_pe:
    input: "mapped_reads/pe/{sample}.phred30.bam"
    output: "mapped_reads/pe/{sample}.sorted.bam"
    conda: "chipseq_env.yaml"
    shell: "samtools sort {input} -o {output}"

# --- POST-MAPPING (PICARD & MACS3) ---

rule remove_duplicates_se:
    input: "mapped_reads/se/{sample}.sorted.bam"
    output: 
        bam="mapped_reads/se/{sample}.noduplicates.bam",
        metrics="mapped_reads/se/{sample}.dup_metrics.txt"
    conda: "chipseq_env.yaml"
    shell: 
        "picard MarkDuplicates I={input} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=STRICT"

rule remove_duplicates_pe:
    input: "mapped_reads/pe/{sample}.sorted.bam"
    output: 
        bam="mapped_reads/pe/{sample}.noduplicates.bam",
        metrics="mapped_reads/pe/{sample}.dup_metrics.txt"
    conda: "chipseq_env.yaml"
    shell: 
        "picard MarkDuplicates I={input} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=STRICT"

rule macs3_se:
    input: bam="mapped_reads/se/{sample}.noduplicates.bam"
    output: "macs3_peaks/se/{sample}_peaks.narrowPeak"
    conda: "chipseq_env.yaml"
    shell: 
        """
        mkdir -p macs3_peaks/se
        macs3 callpeak -t {input.bam} -f BAM -g 2e7 -q 0.001 --nomodel --shift 0 --extsize 200 -n {wildcards.sample} --outdir macs3_peaks/se
        """

rule macs3_pe:
    input: bam="mapped_reads/pe/{sample}.noduplicates.bam"
    output: "macs3_peaks/pe/{sample}_peaks.narrowPeak"
    conda: "chipseq_env.yaml"
    shell: 
        """
        mkdir -p macs3_peaks/pe
        macs3 callpeak -t {input.bam} -f BAM -g 2e7 -q 0.001 --nomodel --shift 0 --extsize 200 -n {wildcards.sample} --outdir macs3_peaks/pe
        """

# --- VISUALIZATION PREP ---

rule chrom_sizes:
    input: genome="ref/P_falciparum3D7.fa"
    output: "chromsizes.genome"
    conda: "chipseq_env.yaml"
    shell:
        """
        samtools faidx {input.genome}
        cut -f1,2 {input.genome}.fai > {output}
        """

rule bed_graph_pe:
    input: "macs3_peaks/pe/{sample}_peaks.narrowPeak"
    output: "bedgraphs/pe/{sample}.bedGraph"
    conda: "chipseq_env.yaml"
    shell: 
        """ 
        mkdir -p bedgraphs/pe
        awk '{{print $1"\t"$2"\t"$3"\t"$7}}' {input} | sort -k1,1 -k2,2n > {output}
        """

rule bed_graph_se:
    input: "macs3_peaks/se/{sample}_peaks.narrowPeak"
    output: "bedgraphs/se/{sample}.bedGraph"
    conda: "chipseq_env.yaml"
    shell: 
        """ 
        mkdir -p bedgraphs/se
        awk '{{print $1"\t"$2"\t"$3"\t"$7}}' {input} | sort -k1,1 -k2,2n > {output}
        """

rule bigwig_pe:
    input:
        sizes="chromsizes.genome",
        bg="bedgraphs/pe/{sample}.bedGraph"
    output: "bigwig_files/pe/{sample}.bw"
    conda: "chipseq_env.yaml"
    shell:
        """
        mkdir -p bigwig_files/pe
        bedGraphToBigWig {input.bg} {input.sizes} {output}
        """

rule bigwig_se:
    input:
        sizes="chromsizes.genome",
        bg="bedgraphs/se/{sample}.bedGraph"
    output: "bigwig_files/se/{sample}.bw"
    conda: "chipseq_env.yaml"
    shell:
        """
        mkdir -p bigwig_files/se
        bedGraphToBigWig {input.bg} {input.sizes} {output}
        """

rule cleanup:
    shell:
        """
        rm -rf ncbi_dataset ncbi_dataset.zip ref data/fastq/se data/fastq/pe trimmed mapped_reads macs3_peaks bedgraphs bigwig_files .snakemake chromsizes.genome
        """
