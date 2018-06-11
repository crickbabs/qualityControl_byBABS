#!/bin/sh

module load GCC/6.4.0-2.28

PROG=qc_directionality

# ----
QCTXTDIR=/camp/stp/babs/working/bahn/code/project/qc_pipeline/txt
RRNA=$QCTXTDIR/human_rrna_ids.txt
MTRRNA=$QCTXTDIR/human_mtrrna_ids.txt
GLOBIN=$QCTXTDIR/human_globin_ids.txt
# ----

# ----
#BED=/camp/stp/babs/working/data/genomes
#BED=$BED/mus_musculus/ensembl/GRCm38/release-86/gtf
#BED=$BED/Mus_musculus.GRCm38.86.bed
# ----

#BAM=/camp/stp/babs/working/bahn/bam/sample01.cutadapt.STAR.genome.sorted.bam

# ----
#BAM=/camp/stp/babs/working/bahn/projects/hayday/maria.iannitto
#BAM=$BAM/009_rnaseq_tigit_thymocyte_maturation/picard/Th38POS.dupmarked.bam
# ----

# ----
BED=/camp/stp/babs/working/data/genomes
BED=$BED/homo_sapiens/ensembl/GRCh38/release-86/gtf
BED=$BED/Homo_sapiens.GRCh38.86.bed
# ----

# ----
#BAM=/camp/stp/babs/working/bahn/code/workflow/rnaseq/nflow
#BAM=$BAM/general/nextflow_single_end/picard/control1.dupmarked.bam
# ----

# ----
#BAM=/camp/stp/babs/working/bahn/code/workflow/rnaseq/nflow
#BAM=$BAM/general/nextflow_paired_end/picard/control1.dupmarked.bam
# ----

# ----
BAM=test/multiqc/paired_end/control1.dupmarked.bam
# ----

# ----
#BAM=test/multiqc/single_end/control1.dupmarked.bam
# ----


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
	./$PROG $BAM $BED $RRNA $MTRRNA $GLOBIN
fi

