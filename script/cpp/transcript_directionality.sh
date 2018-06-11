#!/bin/sh

module load GCC/6.4.0-2.28

PROG=transcript_directionality

## ---
#BED=/camp/stp/babs/working/data/genomes
#BED=/mus_musculus/ensembl/GRCm38/release-86/gtf
#BED=$BED/Mus_musculus.GRCm38.86.bed
#BAM=/camp/stp/babs/working/bahn/bam/sample01.cutadapt.STAR.genome.sorted.bam
## ---

## ---
#BED=/camp/stp/babs/working/data/genomes
#BED=$BED/homo_sapiens/ensembl/GRCh38/release-86/gtf
#BED=$BED/Homo_sapiens.GRCh38.86.bed
#BAM=/camp/stp/babs/working/bahn/projects/hayday/maria.iannitto
#BAM=$BAM/009_rnaseq_tigit_thymocyte_maturation/picard/Th38POS.dupmarked.bam
## ---

## ---
#BED=/camp/stp/babs/working/data/genomes
#BED=$BED/homo_sapiens/ensembl/GRCh38/release-86/gtf
#BED=$BED/Homo_sapiens.GRCh38.86.bed
#BAM=/camp/stp/babs/working/bahn/code/workflow/rnaseq/nflow/general
#BAM=$BAM/nextflow_single_end/picard/control1.dupmarked.bam
## ---

## ---
#BED=/camp/stp/babs/working/data/genomes
#BED=$BED/homo_sapiens/ensembl/GRCh38/release-86/gtf
#BED=$BED/Homo_sapiens.GRCh38.86.bed
#BAM=/camp/stp/babs/working/bahn/code/workflow/rnaseq/nflow/general
#BAM=$BAM/nextflow_paired_end/picard/control1.dupmarked.bam
## ---

## ---
#BED=/camp/stp/babs/working/data/genomes
#BED=$BED/homo_sapiens/ensembl/GRCh38/release-86/gtf
#BED=$BED/Homo_sapiens.GRCh38.86.bed
#BAM=test/multiqc/single_end/control1.dupmarked.bam
## ---

# ---
BED=/camp/stp/babs/working/data/genomes
BED=$BED/homo_sapiens/ensembl/GRCh38/release-86/gtf
BED=$BED/Homo_sapiens.GRCh38.86.bed
BAM=test/multiqc/paired_end/control1.dupmarked.bam
# ---

	#-pthread \
if [[ $1 == "compil" ]]
then
	g++ $PROG.cpp \
		-I /camp/stp/babs/working/bahn/code/cpp/seqan/include \
		-std=c++14 \
		-O3 \
		-DNDEBUG \
		-W \
		-Wall \
		-pedantic \
		-fopenmp \
		-lpthread \
		-lrt \
		-lz \
		-DSEQAN_HAS_ZLIB=1 \
		-DSEQAN_ENABLE_DEBUG=0 \
		-DSEQAN_ENABLE_TESTING=0 \
		-o $PROG
elif [[ $1 == "run" ]]
then
	./$PROG $BAM $BED
fi

