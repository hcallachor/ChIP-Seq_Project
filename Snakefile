#main snakefile
import os
import yaml

#load the yaml file containing the metadata for the samples (paired vs single end)
configfile: "config.yaml"
sample_types = config["samples"]

#samples 
samples = [
    d for d in os.listdir("initial_data")
    if os.path.isdir(os.path.join("initial_data", d))
]

rule all:
    input:
        "ref/P_falciparum3D7.fa",    
        expand("macs2_peaks/{sample}_peaks.narrowPeak", sample=samples),
        expand("results/{sample}.jaccard.txt", sample=samples)

#get paired end fastq files from the sra accessions  
rule fasterq_dump:    
    input:
        "initial_data/{sample}/{sample}.sra"
    output:
        "data/fastq/{sample}_1.fastq",
        "data/fastq/{sample}_2.fastq"
    wildcard_constraints:
        sample="|".join(paired_samples)
    shell:
        """
        mkdir -p data/fastq
        fasterq-dump {input} -O data/fastq --split-files
        """
#get single end fastq files from sra accessions
rule fasterq_dump_single:
    input:
        "initial_data/{sample}/{sample}.sra"
    output:
        "data/fastq/{sample}.fastq"
    wildcard_constraints:
        sample="|".join(single_samples)
    shell:
        """
        mkdir -p data/fastq
        fasterq-dump {input} -O data/fastq
        """

#download the genome and annotations
rule download_reference_genome:
    output:
        genome="ref/P_falciparum3D7.fa",
        gff="ref/P_falciparum3D7_annotations.gff3"
    shell:
        """
        mkdir -p ref
        datasets download genome accession GCF_000002765.6 --include gff3,genome --filename ref/ncbi_dataset.zip
        unzip -o ref/ncbi_dataset.zip -d ref/ncbi_dataset
        cp ref/ncbi_dataset/ncbi_dataset/data/*/*genomic.fna {output.genome}
        cp ref/ncbi_dataset/ncbi_dataset/data/*/*genomic.gff {output.gff}
        rm -rf ref/ncbi_dataset ref/ncbi_dataset.zip
        """

#filter single end reads by quality using Trimmomatic 
rule trimmomatic_se:
    input:
        "data/fastq/{sample}.fastq"
    output:
        "trimmed/{sample}.fastq"
    wildcard_constraints:
        sample="|".join(single_samples)
    shell:
        """ 
        mkdir -p trimmed 
        java -jar ~/chipseq_project/Trimmomatic/trimmomatic-0.39.jar SE -phred33 {input} {output} SLIDINGWINDOW:4:30 MINLEN:35
        """

#filter paired end reads by quality using Trimmomatic 
rule trimmomatic_pe:
    input:
        r1="data/fastq/{sample}_1.fastq",
        r2="data/fastq/{sample}_2.fastq"
    output:
        r1_paired="trimmed/{sample}_1_paired.fastq",
        r1_unpaired="trimmed/{sample}_1_unpaired.fastq",
        r2_paired="trimmed/{sample}_2_paired.fastq",
        r2_unpaired="trimmed/{sample}_2_unpaired.fastq"
    wildcard_constraints:
        sample="|".join(paired_samples)
    shell:
        """
        mkdir -p trimmed 
        java -jar ~/chipseq_project/Trimmomatic/trimmomatic-0.39.jar PE -phred33 {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} SLIDINGWINDOW:4:30 MINLEN:35
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
        fastq="trimmed/{sample}.fastq",
        donecheck="ref/index.done"
    output:
        "mapped_reads/{sample}.bam"
    wildcard_constraints:
        sample="|".join(single_samples)
    shell:
        """
        mkdir -p mapped_reads 
        bwa mem -M {input.genome} {input.fastq} | samtools view -bS - > {output}
        """

#map reads to reference genome using BWA-MEM for paired end reads
rule bwa_mapping_pe:
    input:
        genome="ref/P_falciparum3D7.fa",
        r1="trimmed/{sample}_1_paired.fastq",
        r2="trimmed/{sample}_2_paired.fastq",
        idx="ref/index.done"
    output:
        "mapped_reads/{sample}.bam"
    wildcard_constraints:
        sample="|".join(paired_samples)
    shell:
        """
        mkdir -p mapped_reads
        bwa mem -M {input.genome} {input.r1} {input.r2} | samtools view -bS - > {output}
        """

#quality filter bam files by Phred quality score 30 using samtools
rule filter_by_phred_score:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "mapped_reads/{sample}.phred30.bam"
    shell:
        """
        samtools view -b -q 30 {input} > {output}
        """

#sort bam files by coordinate order using samtools. Picards' MarkDuplicates requires sorted bam files as input
rule sort_bam:
    input:
        "mapped_reads/{sample}.phred30.bam"
    output:
        "mapped_reads/{sample}.sorted.bam"
    shell:
        """
        samtools sort {input} -o {output}
        """

#remove duplicates using Picard’s MarkDuplicates
rule remove_duplicates:
    input:
       "mapped_reads/{sample}.sorted.bam"
    output:
        "mapped_reads/{sample}.noduplicates.bam"
    shell:
        """
        java -jar picard.jar MarkDuplicates I={input} O={output} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=STRICT M=mapped_reads/{wildcards.sample}.dup_metrics.txt
        """

#call peaks using MacS2. Need to input a control I think but idk what it is. idk if the command is 100% correct
rule macs2:
    input:
        bam="mapped_reads/{sample}.noduplicates.bam"
    output:
        "macs2_peaks/{sample}_peaks.narrowPeak"
    shell:
        """
        mkdir -p macs2_peaks
        macs2 callpeak -t {input.bam} -f BAM -g 2e7 -q 0.001 --nomodel --shift 0 --extsize 200 -n {wildcards.sample} --outdir macs2_peaks
        """

#filter out overlapping peaks using bedtools intersect to avoid overcounting
rule bedtools_intersect:
    input:

    output:

    shell:
        """
        bedtools intersect -a <> \ -b <>
        """
#overlay the BED files containing our BED output onto the BED files containing the paper-provided BED output to see where they intersect with pybedtools jaccard
rule pybedtools_jaccard:
    input:
        provided_BED="provided_BED/{sample}.bed",
        our_BED="macs2_peaks/{sample}_peaks.narrowPeak"
    output:
        "results.jaccard"
    shell:
        """
        bedtools jaccard -a <{input.provided_BED[0]}> -b <{input.our_BED[0]}> > results.jaccard
        bedtools jaccard -a <{input.provided_BED[1]}> -b <{input.our_BED[1]}> >> results.jaccard
        bedtools jaccard -a <{input.provided_BED[2]}> -b <{input.our_BED[2]}> >> results.jaccard
        """

#cleanup rule to remove files and run snakemake again
rule cleanup:
    shell:
        """
        rm -rf ncbi_dataset ncbi_dataset.zip
        rm -rf ref
        rm -rf data/fastq
        rm -rf trimmed mapped_reads macs2_peaks
        rm -rf .snakemake
        """
