#main snakefile
import os
import yaml


configfile: 'CompProjectconfig.yaml' #snakemake will pull info from this file

#load the yaml file containing the metadata for the samples (paired vs single end)
single_samples = config["Single End"] #this algins with the formatting of the YAML for SINGLE END READS
paired_samples = config["Paired End"] #this aligns with teh formatting of the YAML for PARIED END READS
samples = config["SRAs"] #all of the samples are going to come in the form of an SRA
reference_genome = config["Reference Genome"]

rule all:
    input:   
        expand("bigwig_files/pe/{sample}.bw", sample=paired_samples),
        expand("bigwig_files/se/{sample}.bw", sample=single_samples)
        

#get paired end fastq files from the sra accessions  
rule fasterq_dump:    
    input:
        "initial_data/{sample}/{sample}.sra" #this is the formatting of the folder made by SampleDownloadPFal.py 
    output:
        "data/fastq/pe/{sample}_1.fastq", # -1 signifies the forward reads
        "data/fastq/pe/{sample}_2.fastq" # -2 signifies the reverse reads 
    shell:
        """
        mkdir -p data/fastq/pe
        fasterq-dump {input} -O data/fastq/pe --split-files
        """
        
#get single end fastq files from sra accessions
rule fasterq_dump_single:
    input:
        "initial_data/{sample}/{sample}.sra" #this is the formatting of the folder made by SampleDownloadPFal.py 
    output:
        "data/fastq/se/{sample}.fastq"
    shell:
        """
        mkdir -p data/fastq/se
        fasterq-dump {input} -O data/fastq/se
        """

'''The sample fasterq-dump codes are designed to make file formatting similar between the single end and paired end'''

#download the genome and annotations
rule download_reference_genome:
    output:
        genome="ref/P_falciparum3D7.fa",
        gff="ref/P_falciparum3D7_annotations.gff3"
    shell:
        """
        mkdir -p ref #making the directory 
        datasets download genome accession {reference_genome} --include gff3,genome --filename ref/ncbi_dataset.zip #ncbi datasets used to get accession
        unzip -o ref/ncbi_dataset.zip -d ref/ncbi_dataset
        cp ref/ncbi_dataset/ncbi_dataset/data/*/*genomic.fna {output.genome}
        cp ref/ncbi_dataset/ncbi_dataset/data/*/*genomic.gff {output.gff}
        rm -rf ref/ncbi_dataset ref/ncbi_dataset.zip
        """

#filter single end reads by quality using Trimmomatic 
rule trimmomatic_se:
    input:
        "data/fastq/se/{sample}.fastq"
    output:
        "trimmed/se/{sample}.fastq"
    shell:
        """ 
        mkdir -p trimmed/se 
        java -jar ~/chipseq_project/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 {input} {output} SLIDINGWINDOW:4:30 MINLEN:35
        """

#filter paired end reads by quality using Trimmomatic 
rule trimmomatic_pe:
    input:
        r1="data/fastq/pe/{sample}_1.fastq",
        r2="data/fastq/pe/{sample}_2.fastq"
    output:
        r1_paired="trimmed/pe/{sample}_1_paired.fastq",
        r1_unpaired="trimmed/pe/{sample}_1_unpaired.fastq",
        r2_paired="trimmed/pe/{sample}_2_paired.fastq",
        r2_unpaired="trimmed/pe/{sample}_2_unpaired.fastq"
    shell:
        """
        mkdir -p trimmed/pe
        java -jar ~/chipseq_project/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} SLIDINGWINDOW:4:30 MINLEN:35
        """

#index reference genome for BWA mapping
rule index_ref:
    input:
        genome="ref/P_falciparum3D7.fa"
    output:
        touch("ref/index.done")
    shell:
        """
        bwa index {input.genome}
        touch {output}
        """

#map reads to reference genome using BWA-MEM for single end reads
rule bwa_mapping_se: #use -M in command to make it compatible with Picard, and pipe command to samtools
    input:
        genome="ref/P_falciparum3D7.fa",
        fastq="trimmed/se/{sample}.fastq",
        donecheck="ref/index.done"
    output:
        "mapped_reads/se/{sample}.bam"
    shell:
        """
        mkdir -p mapped_reads/se 
        bwa mem -M {input.genome} {input.fastq} | samtools view -bS - > {output}
        """

#map reads to reference genome using BWA-MEM for paired end reads
rule bwa_mapping_pe:
    input:
        genome="ref/P_falciparum3D7.fa",
        r1="trimmed/pe/{sample}_1_paired.fastq",
        r2="trimmed/pe/{sample}_2_paired.fastq",
        idx="ref/index.done"
    output:
        "mapped_reads/pe/{sample}.bam"
    shell:
        """
        mkdir -p mapped_reads/pe
        bwa mem -M {input.genome} {input.r1} {input.r2} | samtools view -bS - > {output}
        """

#quality filter SE bam files by Phred quality score 30 using samtools
rule filter_by_phred_score_se:
    input:
        "mapped_reads/se/{sample}.bam"
    output:
        "mapped_reads/se/{sample}.phred30.bam"
    shell:
        """
        samtools view -b -q 30 {input} > {output}
        """

#quality filter PE bam files by Phred quality score 30 using samtools
rule filter_by_phred_score_pe:
    input:
        "mapped_reads/pe/{sample}.bam"
    output:
        "mapped_reads/pe/{sample}.phred30.bam"
    shell:
        """
        samtools view -b -q 30 {input} > {output}
        """

#sort SE bam files by coordinate order using samtools. Picards' MarkDuplicates requires sorted bam files as input
rule sort_bam_se:
    input:
        "mapped_reads/se/{sample}.phred30.bam"
    output:
        "mapped_reads/se/{sample}.sorted.bam"
    shell:
        """
        samtools sort {input} -o {output}
        """

#sort PE bam files by coordinate order using samtools. Picards' MarkDuplicates requires sorted bam files as input
rule sort_bam_pe:
    input:
        "mapped_reads/pe/{sample}.phred30.bam"
    output:
        "mapped_reads/pe/{sample}.sorted.bam"
    shell:
        """
        samtools sort {input} -o {output}
        """

#remove duplicates using Picard’s MarkDuplicates for SE
rule remove_duplicates_se:
    input: 
        "mapped_reads/se/{sample}.sorted.bam"
    output: 
        "mapped_reads/se/{sample}.noduplicates.bam"
    shell: 
        """
        java -jar picard.jar MarkDuplicates I={input} O={output} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=STRICT M=mapped_reads/se/{wildcards.sample}.dup_metrics.txt
        """

#remove duplicates using Picard’s MarkDuplicates for PE
rule remove_duplicates_pe:
    input: 
        "mapped_reads/pe/{sample}.sorted.bam"
    output: 
        "mapped_reads/pe/{sample}.noduplicates.bam"
    shell: 
        """
        java -jar picard.jar MarkDuplicates I={input} O={output} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=STRICT M=mapped_reads/pe/{wildcards.sample}.dup_metrics.txt
        """

#call peaks using MacS3 for SE
rule macs3_se:
    input: 
        bam="mapped_reads/se/{sample}.noduplicates.bam"
    output: 
        "macs3_peaks/se/{sample}_peaks.narrowPeak"
    shell: 
        """
        mkdir -p macs3_peaks/se
        macs3 callpeak -t {input.bam} -f BAM -g 2e7 -q 0.001 --nomodel --shift 0 --extsize 200 -n {wildcards.sample} --outdir macs3_peaks/se
        """

#call peaks using MacS3 for PE
rule macs3_pe:
    input: 
        bam="mapped_reads/pe/{sample}.noduplicates.bam"
    output: 
        "macs3_peaks/pe/{sample}_peaks.narrowPeak"
    shell: 
        """
        mkdir -p macs3_peaks/pe
        macs3 callpeak -t {input.bam} -f BAM -g 2e7 -q 0.001 --nomodel --shift 0 --extsize 200 -n {wildcards.sample} --outdir macs3_peaks/pe
        """

#create the chromosome length file to be used in bigwig file creation
rule chrom_sizes:
    input:
        genome="ref/P_falciparum3D7.fa"
    output:
        "chromsizes.genome"
    shell:
        """
        samtools faidx {input.genome}
        cut -f1,2 {input.genome}.fai > {output}
        """

#convert macs3 narrowpeak output into bedgraph files, to be converted into bigwig files
rule bed_graph_pe:
    input:
        "macs3_peaks/pe/{sample}_peaks.narrowPeak"
    output:
        "bedgraphs/pe/{sample}.bedGraph"
    shell:
        """ 
        mkdir -p bedgraphs/pe
        awk '{{print $1"\t"$2"\t"$3"\t"$7}}' {input} | sort -k1,1 -k2,2n > {sample}.bedGraph
        """

#convert macs3 narrowpeak output into bedgraph files, to be converted into bigwig files
rule bed_graph_se:
    input:
        "macs3_peaks/se/{sample}_peaks.narrowPeak"
    output:
        "bedgraphs/se/{sample}.bedGraph"
    shell:
        """ 
        mkdir -p bedgraphs/se
        awk '{{print $1"\t"$2"\t"$3"\t"$7}}' {input} | sort -k1,1 -k2,2n > {sample}.bedGraph
        """

rule bigwig_pe:
    input:
        "chromsizes.genome",
        "bedgraphs/pe/{sample}.bedGraph"
    output:
        bigwig_files/pe/{sample}.bw
    shell:
        """
        mkdir -p bigwig_files/pe
        bedGraphToBigWig {sample}.bedGraph chromsizes.genome {sample}.bw
        """

rule bigwig_se:
    input:
        "chromsizes.genome",
        "bedgraphs/se/{sample}.bedGraph"
    output:
        bigwig_files/se/{sample}.bw
    shell:
        """
        mkdir -p bigwig_files/se
        bedGraphToBigWig {sample}.bedGraph chromsizes.genome {sample}.bw
        """

#cleanup rule to remove files and run snakemake again
rule cleanup:
    shell:
        """
        rm -rf ncbi_dataset ncbi_dataset.zip
        rm -rf ref
        rm -rf data/fastq/se
        rm -rf data/fastq/pe
        rm -rf trimmed 
        rm -rf mapped_reads 
        rm -rf macs3_peaks
        rm -rf bigwig_files
        rm -rf .snakemake
        """
