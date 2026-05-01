
# Flags Used In the PlasmoPeak Pipeline -- May 2026


## rule fasterq_dump:

-p

-O

-- split-files

## rule fasterq_dump_single:

-p

-O

## rule download_reference_genome:

-p

--include

--filename

-o

-d

-rf

## rule trimmomatic_se:

-p

-phred33

## rule trimmomatic_pe:

-p

-phred33

## rule bwa_mapping_se:

-p

-M

-R

-bS

## rule bwa_mapping_pe:

-p

-M

-R

-bS

## rule filter_by_phred_score_se:

-b

-q

## rule filter_by_phred_score_pe:

-b

-q

## sort_bam_se:

-o

## sort_bam_pe

-o

## rule macs3_se:

-p

-f

-t

-g

-q

--nomodel

--shift

--extsize

-n

--outdir

## rule macs3_pe:

-p

-t

-f

-g

-q

--nomodel

--shift

--extsize

-n

--outdir

## rule chrom_sizes:

-f1,2

## rule bed_graph_pe:

-p

-k1,1

-k2,2n


## rule bed_graph_se:

-p

-k1,1

-k2,2n

## rule bigwig_pe:

-p

## rule bigwig_se:

-p

## rule cleanup

-rf


