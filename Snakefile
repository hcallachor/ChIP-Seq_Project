#main snakefile
import os

#samples (idk if the samples are single or paired end we would have to restructure this if we have both so then trimmomatic will work for each kind of sample)
samples = [
    d for d in os.listdir("initial_data")
    if os.path.isdir(os.path.join("initial_data", d))
]

rule all:
    input:
        expand("macs2_peaks/{sample}_peaks.narrowPeak", sample=samples)

#get fastq files from the sra accessions (need to fix so it does both single and paired end files)
rule fasterq_dump:
    input:
        "initial_data/{sample}/{sample}.sra"
    output:
        "data/fastq/{sample}.fastq"
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

#filter reads by quality using Trimmomatic (need to know if data is paired or single end. SE or PE? if we have both types of samples we need to make another rule that does PE and inputs r1 and r2)
rule trimmomatic:
    input:
        "data/fastq/{sample}.fastq"
    output:
        "trimmed/{sample}.fastq"
    shell:
        """
        mkdir -p trimmed
        trimmomatic SE -phred33 {input} {output} SLIDINGWINDOW:4:30 MINLEN:35
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

#map reads to reference genome using BWA-MEM
rule BWA_mapping: #use -M in command to make it compatible with Picard, and pipe command to samtools
    input:
        genome="ref/P_falciparum3D7.fa",
        fastq="trimmed/{sample}.fastq",
        donecheck="ref/index.done"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        """
        mkdir -p mapped_reads
        bwa mem -M {input.genome} {input.fastq} | samtools view -bS - > {output}
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

#bam files must be binarized prior to chromatin state modeling
#paper specifies bin size 200 which is the default size
rule ChromHMM_binarize: #not sure if I need to specify memory with -mx4000M but google said I need it
    input:
        chromosomelengthfile="chromsizes.genome",
        inputbamdir="",
        cellmarkfiletable=""
    output:
        outputbinarydir=""
    shell:
        """
        java -mx4000M -jar ChromHMM.jar BinarizeBam {input.chromosomelengthfile} {input.inputbamdir} {input.cellmarkfiletable} {output.outputbinardir}
        """
#create the chromosome length file to be used in BinarizeBam
#not 100% sure I did this correct
input:
    genome="ref/P_falciparum3D7.fa"
output:
    "chromsizes.genome"
rule chromlength:
    shell:
        samtools faidx {input.genome} | cut -f1,2 P_falciparum3D7.fa.fai > chromsizes.genome

#map the location of chromatin regulatory states across the genome using ChromHMM 
rule ChromHMM_model:
    input:
        inputdir="",
        assembly=""
    output:
        outputdir=""
    shell:
        """
        java -mx4000M -jar ChromHMM.jar LearnModel {input.inputdir} {output.outputdir} 11 {input.assembly}
        """

#overlay the BED files containing the filtered transcription factor binding sites onto the BED files containing the ChromHMM chromatin state locations to see where they intersect with pybedtools
rule pybedtools:
    input:

    output:

    shell:
        """

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
