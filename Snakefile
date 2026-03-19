#main snakefile

#samples (idk 100% what samples we have and if they are single or paired end we would have to restructure this if we have both so then trimmomatic will work for each kind of sample)
samples = [
"H3K9ac_10h",
"H3K9ac_30h",
"H3K9ac_40h",
"H3K4me1_10h",
"H3K4me1_30h",
"H3K4me1_40h",
"H3K4me3_10h",
"H3K4me3_30h",
"H3K4me3_40h",
"H2AZ_10h",
"H2AZ_30h",
"H2AZ_40h",
"H3K27ac_10h",
"H3K27ac_30h",
"H3K27ac_40h",
"H3K18ac_10h",
"H3K18ac_30h",
"H3K18ac_40h",
"H3K9me3_20h",
"H3K9me3_40h",
"HP1_10h",
"HP1_30h",
"HP1_40h",
"ATAC_10h",
"ATAC_30h",
"ATAC_40h"    
]

rule all:
    input:
        expand("mapped_reads/{sample}.phred30.bam", sample=samples)

#download the genome and annotations
rule download_reference_genome:
    output:
        genome="ref/P_falciparum3D7.fa",
        cds="ref/P_falciparum3D7_cds.fa"
    shell:
        """
        mkdir -p ref
        datasets download genome accession GCF_000002765.6 --include gff3,genome --filename ref/ncbi_dataset.zip
        unzip -o ref/ncbi_dataset.zip -d ref/ncbi_dataset
        cp ref/ncbi_dataset/data/GCF_000002765.6/cds_from_genomic.fna {output.cds}
        cp ref/ncbi_dataset/data/GCF_000002765.6/GCF_000002765.6.fna {output.genome}
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

#call peaks using MacS2
rule MacS2:
    input:
    
    output:

    shell:
        """
        
        """

#map the location of chromatin regulatory states across the genome using ChromHMM 
rule ChromHMM:
    input:

    output:

    shell:
        """

        """

#overlay the BED files containing the filtered transcription factor binding sites onto the BED files containing the ChromHMM chromatin state locations to see where they intersect with pybedtools
rule pybedtools:
    input:

    output:

    shell:
        """
     
        """
