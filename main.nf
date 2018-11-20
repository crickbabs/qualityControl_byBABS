#!/usr/bin/env nextflow

/*
 * harshil.patel@crick.ac.uk
 * nourdine.bah@crick.ac.uk
 * philip.east@crick.ac.uk
 */

// groovy modules
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

// the java modules
import java.nio.file.Paths


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                            FUNCTIONS                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

/* ------------------ */
def is_single_end(dir) {

	// fastq files
	def files = []
	new File(dir).eachFile() { files << it }

	// one or two tags in the filenames ?
	def paths = files.collect{ it.toString() }.sort()

	def tags = paths.collect{
		it.replaceAll(".*_(R\\d)_\\d{3}.fastq.gz", "\$1")
		}.unique().sort()

	single_end = tags.size() == 2 ? false : true

	return single_end
}

/* ---------------------------------- */
def get_rough_read_length(read_length) {
	
	// three babs star indices have been built with these read length parameters
	def starIndexReadLengths = [50, 75, 100]
	
	// take the index with the closest read length to the experiment's
	def diffs = []
	starIndexReadLengths.each() { length ->
		diff = (length - read_length.toInteger()).abs()
		diffs.add(diff)
	}
	def index = diffs.findIndexValues() { i -> i == diffs.min() }[0]
	def rough_read_length = starIndexReadLengths[index.toInteger()]

	return rough_read_length
}

/* ----------------------------------------------------------- */
def get_path(species, filetype, genome_dirpath=GENOME_DIRPATH) {

	def binomial = species.toLowerCase().tokenize(" ")

	if ( binomial[0] == "homo" & binomial[1] == "sapiens" ) {
		annot = HUMAN_ANNOT
	} else if ( binomial[0] == "mus" & binomial[1] == "musculus" ) {
		annot = MOUSE_ANNOT
	}

	def filename = ""
	def directory = ""

	if ( filetype.toLowerCase() == "gtf" ) {

		directory =
			Paths.get(
				genome_dirpath,
				binomial.join("_"),
				"ensembl",
				annot["version"],
				"release-" + annot["release"],
				"gtf"
			).toString()
		
		filename =
			binomial[0].capitalize() +
			"_" +
			binomial[1] +
			"." +
			annot["version"] +
			"." +
			annot["release"] +
			".gtf"
		
	} else if ( filetype.toLowerCase() == "fasta" ) {
		
		directory =
			Paths.get(
				genome_dirpath,
				binomial.join("_"),
				"ensembl",
				annot["version"],
				"release-" + annot["release"],
				"genome"
			).toString()
		
		filename =
			binomial[0].capitalize() +
			"_" +
			binomial[1] +
			"." +
			annot["version"] +
			"." +
			"dna_sm" +
			"." +
			annot["fasta_suffix"] +
			".fa"
		
	} else if ( filetype.toLowerCase() == "bed" ) {
		
		directory =
			Paths.get(
				genome_dirpath,
				binomial.join("_"),
				"ensembl",
				annot["version"],
				"release-" + annot["release"],
				"gtf"
			).toString()
		
		filename =
			binomial[0].capitalize() +
			"_" +
			binomial[1] +
			"." +
			annot["version"] +
			"." +
			annot["release"] +
			".bed"
		
	} else if ( filetype.toLowerCase() == "rnaseqc_gtf" ) {
		
		directory =
			Paths.get(
				genome_dirpath,
				binomial.join("_"),
				"ensembl",
				annot["version"],
				"release-" + annot["release"],
				"gtf"
			).toString()
		
		filename =
			binomial[0].capitalize() +
			"_" +
			binomial[1] +
			"." +
			annot["version"] +
			"." +
			annot["release"] +
			".rnaseqc.gtf"
		
	} else if ( filetype.toLowerCase() == "refflat" ) {
		
		directory =
			Paths.get(
				genome_dirpath,
				binomial.join("_"),
				"ensembl",
				annot["version"],
				"release-" + annot["release"],
				"gtf"
			).toString()
		
		filename =
			binomial[0].capitalize() +
			"_" +
			binomial[1] +
			"." +
			annot["version"] +
			"." +
			annot["release"] +
			".refflat"
		
	} else if ( filetype.toLowerCase() == "rrna" ) {
		
		directory =
			Paths.get(
				genome_dirpath,
				binomial.join("_"),
				"ensembl",
				annot["version"],
				"release-" + annot["release"],
				"gtf"
			).toString()
		
		filename =
			binomial[0].capitalize() +
			"_" +
			binomial[1] +
			"." +
			annot["version"] +
			"." +
			annot["release"] +
			".rRNA.list"
		
	} else if ( filetype.toLowerCase() == "rrna_interval" ) {
		
		directory =
			Paths.get(
				genome_dirpath,
				binomial.join("_"),
				"ensembl",
				annot["version"],
				"release-" + annot["release"],
				"gtf"
			).toString()
		
		filename =
			binomial[0].capitalize() +
			"_" +
			binomial[1] +
			"." +
			annot["version"] +
			"." +
			annot["release"] +
			".rRNA.interval_list"
	}

	return Paths.get(directory, filename).toString()
}

/* ----------------------------- */
def get_idx_path(species, length) {

	def binomial = species.toLowerCase().tokenize(" ")

	if ( binomial[0] == "homo" & binomial[1] == "sapiens" ) {
		annot = HUMAN_ANNOT
	} else if ( binomial[0] == "mus" & binomial[1] == "musculus" ) {
		annot = MOUSE_ANNOT
	}

	def path =
		Paths.get(
			GENOME_DIRPATH,
			binomial.join("_"),
			"ensembl",
			annot["version"],
			"release-" + annot["release"],
			"genome_idx",
			"rsem",
			"star",
			length + "bp",
			"genome"
		).toString()

	return path
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         SCRIPTS PATHS                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////


def CHROM_SCRIPT =
	Paths.get(
		workflow.projectDir.toString(),
		"scripts",
		"py",
		"reads_per_chrom.py"
		)

def QC_GENES_SCRIPT =
	Paths.get(
		workflow.projectDir.toString(),
		"scripts",
		"cpp",
		"qc_directionality"
		)

def TRANSCRIPTS_SCRIPT =
	Paths.get(
		workflow.projectDir.toString(),
		"scripts",
		"cpp",
		"transcript_directionality"
		)
def TXT_DIRPATH =
	Paths.get(
		workflow.projectDir.toString(),
		"scripts",
		"qc_genes_list",
		"txt"
		).toString()
def RRNA = Paths.get(TXT_DIRPATH, "human_rrna_ids.txt")
def MTRRNA = Paths.get(TXT_DIRPATH, "human_mtrrna_ids.txt")
def GLOBIN = Paths.get(TXT_DIRPATH, "human_globin_ids.txt")

def TSS_SCRIPT =
	Paths.get(
		workflow.projectDir.toString(),
		"scripts",
		"py",
		"tss_coverage.py"
		)


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                              MODULES                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

def MODULE_ANACONDA = "Anaconda2/5.1.0"
def MODULE_CUTADAPT = "cutadapt/1.9.1-foss-2016b-Python-2.7.12"
def MODULE_FASTQC = "FastQC/0.11.7-Java-1.8.0_162"
def MODULE_FSCREEN = "fastq_screen/0.9.3-2016a-Perl-5.22.1"
def MODULE_GCC = "GCC/6.4.0-2.28"
def MODULE_HTSEQ = "HTSeq/0.6.1p1-foss-2016b-Python-2.7.12"
def MODULE_MULTIQC = "multiqc/1.3-2016b-Python-2.7.12"
def MODULE_PANDAS = "pandas/0.16.2-2016b-Python-2.7.12"
def MODULE_PICARD = "picard/2.1.1-Java-1.8.0_92"
def MODULE_PYTHON = "Python/3.6.3-foss-2017b"
def MODULE_R = "R/3.4.0-intel-2017a-X11-20170314"
def MODULE_RNASEQC = "RNA-SeQC/1.1.8-Java-1.7.0_80"
def MODULE_RSEM = "RSEM/1.3.0-foss-2016b"
def MODULE_RSEQC = "RSeQC/2.6.4-foss-2016b-Python-2.7.12-R-3.3.1"
def MODULE_SAMTOOLS = "SAMtools/1.3.1-foss-2016b"
def MODULE_SEQTK = "seqtk/1.2-foss-2016b"
def MODULE_STAR = "STAR/2.5.2a-foss-2016b"


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       CONFIGURATION FILES                           -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

def FSCREEN_CONF_FILEPATH =
	Paths.get(workflow.projectDir.toString(),
	"conf",
	"fastq_screen.conf"
	).toString()

def MULTIQC_CONF_FILEPATH =
	Paths.get(workflow.projectDir.toString(),
	"multiqc",
	"conf.yml"
	).toString()


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       OUTPUT DIRECTORIES                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
 
def OUTPUT_DIRNAME = "output"

def CHROM_DIRNAME = "chromosome"
def CUTADAPT_DIRNAME = "cutadapt"
def FSCREEN_DIRNAME = "fscreen"
def IDXSTATS_DIRNAME = "idxstats"
def INFO_DIRNAME = "info"
def QC_GENES_DIRNAME = "qc_genes"
def RNASEQC_DIRNAME = "rnaseqc"
def SAMPLING_DIRNAME = "sampling"
def STAR_DIRNAME = "rsem_star"
def STATS_DIRNAME = "stats"
def TRANSCRIPTS_DIRNAME = "transcripts"
def TSS_DIRNAME = "tss"

def FASTQC_DIRNAME = "fastqc"
def FASTQC_RAW_DIRNAME = "raw"
def FASTQC_CUTADAPT_DIRNAME = "cutadapt"

def PICARD_DIRNAME = "picard"
def PICARD_GROUP_DIRNAME = "group"
def PICARD_DUPLICATE_DIRNAME = "duplicate"
def PICARD_COMPLEXITY_DIRNAME = "complexity"
def PICARD_RNASEQMETRICS_DIRNAME = "rnaseqmetrics"
def PICARD_MULTIMETRICS_DIRNAME = "multimetrics"

def RSEQC_DIRNAME = "rseqc"
def RSEQC_INFER_EXPERIMENT_DIRNAME = "infer_experiment"
def RSEQC_JUNCTION_ANNOTATION_DIRNAME = "junction_annotation"
def RSEQC_JUNCTION_SATURATION_DIRNAME = "junction_saturation"
def RSEQC_MISMATCH_PROFILE_DIRNAME = "mismatch_profile"
def RSEQC_READ_DISTRIBUTION_DIRNAME = "read_distribution"
def RSEQC_TRANSCRIPT_INTEGRITY_DIRNAME = "transcript_integrity"

def OUTPUT_DIRPATH =
	Paths.get(workflow.projectDir.toString(), OUTPUT_DIRNAME).toString()

def CHROM_DIRPATH = Paths.get(OUTPUT_DIRPATH, CHROM_DIRNAME).toString()
def CUTADAPT_DIRPATH = Paths.get(OUTPUT_DIRPATH, CUTADAPT_DIRNAME).toString()
def FSCREEN_DIRPATH = Paths.get(OUTPUT_DIRPATH, FSCREEN_DIRNAME).toString()
def IDXSTATS_DIRPATH = Paths.get(OUTPUT_DIRPATH, IDXSTATS_DIRNAME).toString()
def INFO_DIRPATH = Paths.get(OUTPUT_DIRPATH, INFO_DIRNAME).toString()
def QC_GENES_DIRPATH = Paths.get(OUTPUT_DIRPATH, QC_GENES_DIRNAME).toString()
def RNASEQC_DIRPATH = Paths.get(OUTPUT_DIRPATH, RNASEQC_DIRNAME).toString()
def SAMPLING_DIRPATH = Paths.get(OUTPUT_DIRPATH, SAMPLING_DIRNAME).toString()
def STAR_DIRPATH = Paths.get(OUTPUT_DIRPATH, STAR_DIRNAME).toString()
def STATS_DIRPATH = Paths.get(OUTPUT_DIRPATH, STATS_DIRNAME).toString()
def TRANSCRIPTS_DIRPATH =
	Paths.get(OUTPUT_DIRPATH, TRANSCRIPTS_DIRNAME).toString()
def TSS_DIRPATH = Paths.get(OUTPUT_DIRPATH, TSS_DIRNAME).toString()

def FASTQC_DIRPATH = Paths.get(OUTPUT_DIRPATH, FASTQC_DIRNAME).toString()
def FASTQC_RAW_DIRPATH =
	Paths.get(FASTQC_DIRPATH, FASTQC_RAW_DIRNAME).toString()
def FASTQC_CUTADAPT_DIRPATH =
	Paths.get(FASTQC_DIRPATH, FASTQC_CUTADAPT_DIRNAME).toString()

def PICARD_DIRPATH = Paths.get(OUTPUT_DIRPATH, PICARD_DIRNAME).toString()
def PICARD_GROUP_DIRPATH =
	Paths.get(PICARD_DIRPATH, PICARD_GROUP_DIRNAME).toString()
def PICARD_DUPLICATE_DIRPATH =
	Paths.get(PICARD_DIRPATH, PICARD_DUPLICATE_DIRNAME).toString()
def PICARD_COMPLEXITY_DIRPATH =
	Paths.get(PICARD_DIRPATH, PICARD_COMPLEXITY_DIRNAME).toString()
def PICARD_RNASEQMETRICS_DIRPATH =
	Paths.get(PICARD_DIRPATH, PICARD_RNASEQMETRICS_DIRNAME).toString()
def PICARD_MULTIMETRICS_DIRPATH =
	Paths.get(PICARD_DIRPATH, PICARD_MULTIMETRICS_DIRNAME).toString()

def RSEQC_DIRPATH = Paths.get(OUTPUT_DIRPATH, RSEQC_DIRNAME).toString()
def RSEQC_INFER_EXPERIMENT_DIRPATH =
	Paths.get(RSEQC_DIRPATH, RSEQC_INFER_EXPERIMENT_DIRNAME).toString()
def RSEQC_JUNCTION_ANNOTATION_DIRPATH =
	Paths.get(RSEQC_DIRPATH, RSEQC_JUNCTION_ANNOTATION_DIRNAME).toString()
def RSEQC_JUNCTION_SATURATION_DIRPATH =
	Paths.get(RSEQC_DIRPATH, RSEQC_JUNCTION_SATURATION_DIRNAME).toString()
def RSEQC_MISMATCH_PROFILE_DIRPATH =
	Paths.get(RSEQC_DIRPATH, RSEQC_MISMATCH_PROFILE_DIRNAME).toString()
def RSEQC_READ_DISTRIBUTION_DIRPATH =
	Paths.get(RSEQC_DIRPATH, RSEQC_READ_DISTRIBUTION_DIRNAME).toString()
def RSEQC_TRANSCRIPT_INTEGRITY_DIRPATH =
	Paths.get(RSEQC_DIRPATH, RSEQC_TRANSCRIPT_INTEGRITY_DIRNAME).toString()


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          PARAMETERS                                 -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

def SPECIES = params.species
def STRANDEDNESS = params.strandedness
def DIRECTORY = params.directory


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        VARIOUS SETTINGS                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// be careful, no def because it's a global variable
GENOME_DIRPATH = "/camp/stp/babs/working/data/genomes"

def ANACONDA_ENV = "/camp/stp/babs/working/software/anaconda/envs/qc_pipeline"

def PUBLISHDIR_MODE = "copy"
def PUBLISHDIR_OVERWRITE = true


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     SPECIES ANNOTATIONS                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// be careful, no def because it's a global variable
HUMAN_ANNOT = [:]
HUMAN_ANNOT["version"] = "GRCh38"
HUMAN_ANNOT["release"] = "86"
HUMAN_ANNOT["fasta_suffix"] = "primary_assembly"

// be careful, no def because it's a global variable
MOUSE_ANNOT = [:]
MOUSE_ANNOT["version"] = "GRCm38"
MOUSE_ANNOT["release"] = "86"
MOUSE_ANNOT["fasta_suffix"] = "primary_assembly"


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          IS SINGLE-END ?                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

def SINGLE_END = is_single_end(DIRECTORY)


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       THE INPUT CHANNELS                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// transmit to read length determination and to fastqc
Channel
	.fromPath( DIRECTORY + "/*.fastq.gz" )
	.into{
		read_length_fastq;
		fastqc_fastq
		}

if (SINGLE_END) {

	samples =
		Channel
		.fromPath( DIRECTORY + "/*.fastq.gz" )
		.map{[
			"name":it.toString().replaceAll(".*\\/(.*)\\.fastq.gz", "\$1"),
			"file": it
			]}

} else {

	samples =
		Channel
		.fromFilePairs( DIRECTORY + "/*_R{1,2}_*.fastq.gz" )
		.map{[
			"name":it[0]
						.toString()
						.replaceAll(".*\\/(.*)_R\\d_\\d{3}\\.fastq.gz", "\$1"),
			"file1": it[1][0], "file2": it[1][1]
			]}
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                           READ LENGTH                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// get read length from fastq files
process fastq_read_length {

	// custom label
	tag { fastq }

	// problem with slurm otherwise
	beforeScript "module purge"

	// HPC
	cpus 1
	executor "local"
	memory "1000"

	input:
		file(fastq) from read_length_fastq

	output:
		file "*.read_length" into fastq_read_length

	shell:
		"""
		zcat ${fastq} \
			| head -n 16 \
			| sed -n "2~4p" \
			| awk '{ print length }' \
			| sort \
			| uniq \
			| sed -n "1p" \
			> ${fastq}.read_length
		"""
}

// the mean of read lengths for all fastq files
process read_length {

	// custom label
	tag { lengths }

	// problem with slurm otherwise
	beforeScript "module purge"

	// HPC
	cpus 1
	executor "local"
	memory "1000"

	// output directory
	publishDir INFO_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		file(lengths) from fastq_read_length.collect()

	output:
		file("read_length.txt") into read_length

	shell:
		"""
		N=\$(ls ${lengths} | wc -l)
		cat ${lengths} \
			| tr "\n" "+" \
			| sed 's/+\$//' \
			| awk -v N="\$N" '{ print "scale=0;(" \$0 ")/" N }' \
			| bc \
			> read_length.txt
		"""
}

// get read length from the file
read_length.splitText().set{ read_length }
READ_LENGTH = read_length.collect().get().getAt(0)
ROUGH_READ_LENGTH = get_rough_read_length(READ_LENGTH)


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                      ANNOTATION FILE PATHS                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

def ANNOT_GTF_FILEPATH = get_path(SPECIES, "gtf")
def SEQ_FILEPATH = get_path(SPECIES, "fasta")
def ANNOT_BED_FILEPATH = get_path(SPECIES, "bed")
def ANNOT_RNASEQC_FILEPATH = get_path(SPECIES, "rnaseqc_gtf")
def ANNOT_REFFLAT_FILEPATH = get_path(SPECIES, "refflat")
def ANNOT_RRNA_FILEPATH = get_path(SPECIES, "rrna")
def ANNOT_RRNA_INTERVAL_FILEPATH = get_path(SPECIES, "rrna_interval")

RSEM_STAR_INDICE_PREFIX_FILEPATH = get_idx_path(SPECIES, ROUGH_READ_LENGTH)

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                              FASTQC                                 -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process fastqc {

	// custom label
	tag { name }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_FASTQC

	// HPC
	cpus 1
	executor "slurm"
	memory "6000"

	// output directory
	publishDir FASTQC_RAW_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		file fastq from fastqc_fastq

	output:
		file "*fastqc*" into fastqc
	
	shell:
	
		name = fastq.toString().replaceFirst(".fastq.gz", "")
	
		"""
		fastqc ${fastq}
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             CUTADAPT                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process cutadapt {

	// custom label
	tag { name }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_CUTADAPT

	// HPC
	cpus 1
	executor "slurm"
	memory "6000"

	// output directory
	publishDir CUTADAPT_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		val sample from samples
	
	output:
		set val(name), file('*.log') into cutadapt_log
		set val(name), file('*.cutadapt.fastq.gz') \
			into \
				cutadapt_sampling,
				cutadapt_star,
				cutadapt_fscreen,
				cutadapt_fastqc,
				cutadapt_multiqc

	////////////////////////////////////////////////////////////////////////////
	shell:

		// sample name
		name = sample["name"]

		if (SINGLE_END) {

			// locations
			fastq = sample["file"]

			// command arguments
			adapters = "-a AGATCGGAAGAGC"
			input = "-o " + name + ".cutadapt.fastq.gz"
			output = fastq

		} else {

			// locations
			name_1 = sample["name"] + "_1.cutadapt.fastq.gz"
			name_2 = sample["name"] + "_2.cutadapt.fastq.gz"
			fastq_1 = sample["file1"]
			fastq_2 = sample["file2"]

			// command arguments
			adapters = "-a AGATCGGAAGAGC -A AGATCGGAAGAGC"
			input = "-o " + name_1 + " -p " + name_2
			output = fastq_1 + " " + fastq_2
		}

		"""
		cutadapt \
			${adapters} \
			${input} \
			-e 0.1 \
			-q 10 \
			-m 25 \
			-O 1 \
			${output} > ${name}.log
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         FASTQC CUTADAPT                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process fastqc_cutadapt {

	// custom label
	tag { name }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_FASTQC

	// HPC
	cpus 1
	executor "slurm"
	memory "6000"

	// output directory
	publishDir FASTQC_CUTADAPT_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE


	input:
		set val(name), file(fastq) from cutadapt_fastqc

	output:
		file "*fastqc*" into fastqc_cutadapt
	
	shell:
	
		name = fastq.toString().replaceFirst(".fastq.gz", "")
	
		"""
		fastqc ${fastq}
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --			                     SAMPLING                                 -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process sampling {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_SEQTK

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir SAMPLING_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(fastq) from cutadapt_sampling

	output:
		set val(sample), file("*.fastq.gz") into sampling, sampling_multiqc
	
	shell:

		n = 1000
		
		if (SINGLE_END) {
		
			name = fastq.toString().replaceFirst(".fastq.gz", "")
		
			"""
			seqtk sample ${fastq} ${n} | gzip -c > ${name}.sampled.fastq.gz
			"""
		
		} else {
		
			name = fastq.collect { it.toString().replaceFirst(".fastq.gz", "") }
		
			"""
			seqtk sample -s1903 ${fastq[0]} ${n} \
				| gzip -c > ${name[0]}.sampled.fastq.gz
			seqtk sample -s1903 ${fastq[1]} ${n} \
				| gzip -c > ${name[1]}.sampled.fastq.gz
			"""
		}
}

process sampled_rsem_star {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEM
	module MODULE_STAR

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir SAMPLING_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(fastq) from sampling

	output:
		set val(sample), file("*.STAR.genome.bam") \
			into \
				sampled_bam,
				sampled_bam_multiqc

	shell:

		if (SINGLE_END) {

			input = fastq

		} else {

			input = "--paired-end " + fastq[0] + " " + fastq[1]
		}

		"""
		rsem-calculate-expression \
			--temporary-folder "tmp" \
			--star \
			--num-threads ${task.cpus} \
			--strandedness ${STRANDEDNESS} \
			--estimate-rspd \
			--seed 1 \
			--output-genome-bam \
			--star-output-genome-bam \
			--star-gzipped-read-file \
			${input} \
			${RSEM_STAR_INDICE_PREFIX_FILEPATH} \
			${sample}
		"""
}

process sampled_sort {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_SAMTOOLS

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir SAMPLING_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from sampled_bam

	output:
		set val(sample), file("*.bam"), file("*.bai") \
			into \
				sampled_sorted_bam,
				sampled_sorted_bam_multiqc

	shell:
		filename = sample + ".sorted.bam"
		"""
		samtools sort \
			--threads ${task.cpus} \
			-o ${filename} \
			${bam}
		samtools index ${filename}
		"""
}

process sampled_picard_group {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_PICARD

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir SAMPLING_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam), file(bai) from sampled_sorted_bam

	output:
		set val(sample), file("*.bam") \
			into \
				sampled_picard_rg,
				sampled_picard_rg_multiqc

	shell:

		tmp_dirname = "tmp"
		filename = sample + ".rg.bam"

		"""
		java -Xmx10g -Djava.io.tmpdir=${tmp_dirname} \
			-jar \$EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
			VALIDATION_STRINGENCY=SILENT \
			INPUT=${bam} \
			OUTPUT=${filename} \
			RGID=${sample} \
			RGLB=${sample} \
			RGPU=${sample} \
			RGSM=${sample} \
			RGCN=TheFrancisCrickInsitute \
			RGPL=Illumina
		"""
}

process sampled_picard_duplicate {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_PICARD

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir SAMPLING_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from sampled_picard_rg

	output:
		set val(sample), file("*.bam") \
			into \
				sampled_picard_duplicate,
				sampled_picard_duplicate_multiqc

	shell:

		tmp_dirname = "tmp"
		filename = sample + ".dupmarked.bam"

		// the metrics file
		metrics_filename = sample + ".marked_duplicates"

		"""
		java -Xmx10g -Djava.io.tmpdir=${tmp_dirname} \
			-jar \$EBROOTPICARD/picard.jar MarkDuplicates \
			VALIDATION_STRINGENCY=SILENT \
			INPUT=${bam} \
			OUTPUT=${filename} \
			METRICS_FILE=${metrics_filename} \
			ASSUME_SORTED=true \
			REMOVE_DUPLICATES=false \
			TMP_DIR=${tmp_dirname}
		"""
}

process sampled_picard_index {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_SAMTOOLS

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir SAMPLING_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from sampled_picard_duplicate

	output:
		set val(sample), file(bam), file("*.bai") \
			into\
				sampled_bam_infer_experiment,
				sampled_bam_rnaseqc,
				sampled_picard_index_multiqc

	shell:
		"""
		samtools index ${bam}
		"""
}

process sampled_infer_experiment {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir SAMPLING_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam), file(bai) from sampled_bam_infer_experiment
	
	output:
		set val(sample), file("*.infer_experiment*") \
			into \
				sampled_infer_experiment,
				sampled_infer_experiment_multiqc
	
	shell:

		metrics_filename = sample + ".infer_experiment"

		"""
		infer_experiment.py \
			-i ${bam} \
			-r ${ANNOT_BED_FILEPATH} \
			> ${metrics_filename}
		"""
}

process sampled_rnaseqc {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RNASEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir SAMPLING_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam), file(bai) from sampled_bam_rnaseqc

	output:
		set val(sample), file("*rnaseqc*") \
			into \
				sampled_rnaseqc,
				sampled_rnaseqc_multiqc

	shell:
		
		metrics_filename = sample + ".rnaseqc"

		if (SINGLE_END) {

			single_end = "-singleEnd"

		} else {

			single_end = ""

		}

		"""
		java -Xmx10g -jar \$EBROOTRNAMINSEQC/RNA-SeQC_v[0-9].[0-9].[0-9].jar \
			-d 1000000 \
			-rRNA ${ANNOT_RRNA_FILEPATH} \
			-r ${SEQ_FILEPATH} \
			-t ${ANNOT_RNASEQC_FILEPATH} \
			-o ${metrics_filename} \
			-gatkFlags '-S SILENT -U ALLOW_SEQ_DICT_INCOMPATIBILITY' \
			${single_end} \
			-s "${sample}|${bam}|${sample}"
		"""
}

process sampled_multiqc {

	echo true

	// problem with slurm otherwise
	beforeScript "module purge"

	// module
	//module MODULE_PYTHON
	//module MODULE_MULTIQC
	module MODULE_ANACONDA

	// anaconda
	conda ANACONDA_ENV

	// HPC
	cpus 1
	executor "slurm"
	memory "6000"

	// output directory
	publishDir SAMPLING_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		file("*") from sampled_infer_experiment_multiqc.map{it[1]}.collect()
		file("*") from sampled_rnaseqc_multiqc.map{it[1]}.collect()

	output:
		file "multiqc_data/multiqc_rseqc_infer_experiment.json" into rseqc_json
		file "multiqc_data/multiqc_rna_seqc.json" into rnaseqc_json

	shell:
		"""
		export LC_ALL=en_US.utf8
		export LANG=en_US.utf8
		multiqc \
			--data-format json \
			--template custom \
			--config ${MULTIQC_CONF_FILEPATH} \
			.
		"""
}

process strandedness {

	// HPC
	cpus 1
	executor "slurm"
	memory "6000"

	// output directory
	publishDir SAMPLING_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		file rseqc from rseqc_json
		file rnaseqc from rnaseqc_json
	
	output:
		file "strandedness.txt" \
			into \
				strandedness,
				strandedness_multiqc,
				strandedness_variable
	
	shell:
		
		if (SINGLE_END) {

			"""
			#!/usr/bin/python

			import json
			
			rseqc_filepath = "multiqc_rseqc_infer_experiment.json"
			rnaseqc_filepath = "multiqc_rna_seqc.json"
			
			# load json as dict
			with open(rseqc_filepath, 'r') as f:
				rseqc = json.load(f)
			with open(rnaseqc_filepath, 'r') as f:
				rnaseqc = json.load(f)
			
			# rseqc
			rseqc_results = []
			for k, v in rseqc.items():
				# we are getting the index with the alphabetical order (["failed",
				# "se/pe_antisense", "se/pe_sense"])
				keys = sorted(list( v.keys() ))
				sense = v[ keys[2] ]
				antisense = v[ keys[1] ]
				nrmsd_antisense = antisense / (sense + antisense) * 100.0
				nrmsd_sense = sense / (sense + antisense) * 100.0
				rseqc_results.append(nrmsd_sense)
			rseqc_result = sum(rseqc_results) / len(rseqc_results)
			
			# threshold
			if rseqc_result < 30:
				strandedness = "reverse"
			elif rseqc_result > 80:
				strandedness = "forward"
			else:
				strandedness = "none"

			# result as file
			with open("strandedness.txt", 'w') as f:
				f.write(strandedness)
			"""

		} else {

			"""
			#!/usr/bin/python

			import json
			
			rseqc_filepath = "multiqc_rseqc_infer_experiment.json"
			rnaseqc_filepath = "multiqc_rna_seqc.json"
			
			# load json as dict
			with open(rseqc_filepath, 'r') as f:
				rseqc = json.load(f)
			with open(rnaseqc_filepath, 'r') as f:
				rnaseqc = json.load(f)
			
			# rseqc
			rseqc_results = []
			for k, v in rseqc.items():
				# we are getting the index with the alphabetical order (["failed",
				# "se/pe_antisense", "se/pe_sense"])
				keys = sorted(list( v.keys() ))
				sense = v[ keys[2] ]
				antisense = v[ keys[1] ]
				nrmsd_antisense = antisense / (sense + antisense) * 100.0
				nrmsd_sense = sense / (sense + antisense) * 100.0
				rseqc_results.append(nrmsd_sense)
			rseqc_result = sum(rseqc_results) / len(rseqc_results)
			
			## rnaseqc
			#rnaseqc_results = []
			#for k, v in rnaseqc.items():
			#	res = ( v["End 2 % Sense"] + v["End 2 % Sense"] ) / 2.0
			#	rnaseqc_results.append(res)
			#rnaseqc_result = sum(rnaseqc_results) / len(rnaseqc_results)

			# threshold
			if rseqc_result < 30:
				strandedness = "reverse"
			elif rseqc_result > 80:
				strandedness = "forward"
			else:
				strandedness = "none"

			# result as file
			with open("strandedness.txt", 'w') as f:
				f.write(strandedness)
			"""
		}
}

//// get inferred strandedness from the file
//strandedness_variable.splitText().set{ strandedness_variable }
//INFERRED_STRANDEDNESS = strandedness_variable.collect().get().getAt(0)


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                              INFO                                   -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process info {

	// problem with slurm otherwise
	beforeScript "module purge"

	// HPC
	cpus 1
	executor "slurm"
	memory "1000"

	// output directory
	publishDir INFO_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE
	
	output:
		file "experiment_info.json" into info, info_qc_genes, info_transcripts

	shell:

		def info = [
			protocol: SINGLE_END ? "single-end" : "paired-end",
			strandedness: STRANDEDNESS,
			read_length: READ_LENGTH,
			species: SPECIES
		]

		def json = JsonOutput.toJson(info)

		"""
		echo -E '${json}' > experiment_info.json
		"""
}

process cp_info_transcripts {

	// problem with slurm otherwise
	beforeScript "module purge"

	// HPC
	cpus 1
	executor "slurm"
	memory "1000"

	// output directory
	publishDir TRANSCRIPTS_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE
	
	input:
		file info from info_transcripts
	
	output:
		file "experiment_info_transcripts.json" into info_transcripts_multiqc

	shell:
		"""
		cp -v ${info} experiment_info_transcripts.json
		"""
}

process cp_info_qc_genes {

	// problem with slurm otherwise
	beforeScript "module purge"

	// HPC
	cpus 1
	executor "slurm"
	memory "1000"

	// output directory
	publishDir QC_GENES_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE
	
	input:
		file info from info_qc_genes
	
	output:
		file "experiment_info_qc_genes.json" into info_qc_genes_multiqc

	shell:
		"""
		cp -v ${info} experiment_info_qc_genes.json
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --			                  FASTQ SCREEN                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process fscreen {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_FSCREEN

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir FSCREEN_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(fastq) from cutadapt_fscreen

	output:
		set val(sample), file("*.html") into fscreen_html
		set val(sample), file("*.txt") into fscreen_txt
	
	shell:
		"""
		fastq_screen \
			--force \
			--outdir ./ \
			--subset 200000 \
			--conf ${FSCREEN_CONF_FILEPATH} \
			--threads ${task.cpus} \
			--aligner bowtie2 \
			${fastq}
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                               STAR                                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process rsem_star {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEM
	module MODULE_STAR

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir STAR_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		file(strandedness) from strandedness
		set val(sample), file(fastq) from cutadapt_star

	output:
		set val(sample), file("*.transcript.bam") into rsem_transcript
		set val(sample), file("*.results") into rsem_results
		set val(sample), file("*.stat") into rsem_stat
		set val(sample), file("*.STAR.genome.bam") \
			into \
				rsem_star_genome_sort,
				rsem_star_genome_multiqc

	////////////////////////////////////////////////////////////////////////////
	shell:

		if (SINGLE_END) {

			input = fastq

		} else {

			input = "--paired-end " + fastq[0] + " " + fastq[1]
		}

		"""
		STRANDEDNESS=\$(cat ${strandedness})
		rsem-calculate-expression \
			--temporary-folder "tmp" \
			--star \
			--num-threads ${task.cpus} \
			--strandedness \$STRANDEDNESS \
			--estimate-rspd \
			--seed 1 \
			--output-genome-bam \
			--star-output-genome-bam \
			--star-gzipped-read-file \
			${input} \
			${RSEM_STAR_INDICE_PREFIX_FILEPATH} \
			${sample}
		"""
}

process sort_index_star {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_SAMTOOLS

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir STAR_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from rsem_star_genome_sort

	output:
		set val(sample), file("*.bam"), file("*.bai") \
			into \
				sorted_bam_picard,
				sorted_bam_multiqc

	shell:
		filename = sample + ".sorted.bam"
		"""
		samtools sort \
			--threads ${task.cpus} \
			-o ${filename} \
			${bam}
		samtools index ${filename}
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                            PICARD                                   -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process picard_group {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_PICARD

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir PICARD_GROUP_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam), file(bai) from sorted_bam_picard

	output:
		set val(sample), file("*.bam") into picard_rg_duplicate,
			picard_rg_multiqc

	shell:

		tmp_dirname = "tmp"
		filename = sample + ".rg.bam"

		"""
		java -Xmx10g -Djava.io.tmpdir=${tmp_dirname} \
			-jar \$EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
			VALIDATION_STRINGENCY=SILENT \
			INPUT=${bam} \
			OUTPUT=${filename} \
			RGID=${sample} \
			RGLB=${sample} \
			RGPU=${sample} \
			RGSM=${sample} \
			RGCN=TheFrancisCrickInsitute \
			RGPL=Illumina
		"""
}

process picard_duplicate {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_PICARD

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir PICARD_DUPLICATE_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_rg_duplicate

	output:
		set val(sample), file("*.marked_duplicates") into picard_duplicate
		set val(sample), file("*.bam") \
			into \
				picard_bam_index,
				picard_bam_complexity,
				picard_bam_rnaseqmetrics,
				picard_bam_multimetrics,
				picard_bam_infer_experiment,
				picard_bam_junction_annotation,
				picard_bam_junction_saturation,
				picard_bam_mismatch_profile,
				picard_bam_read_distribution,
				picard_bam_rnaseqc,
				picard_bam_qc_genes,
				picard_bam_transcripts

	shell:

		tmp_dirname = "tmp"
		filename = sample + ".dupmarked.bam"

		// the metrics file
		metrics_filename = sample + ".marked_duplicates"

		"""
		java -Xmx10g -Djava.io.tmpdir=${tmp_dirname} \
			-jar \$EBROOTPICARD/picard.jar MarkDuplicates \
			VALIDATION_STRINGENCY=SILENT \
			INPUT=${bam} \
			OUTPUT=${filename} \
			METRICS_FILE=${metrics_filename} \
			ASSUME_SORTED=true \
			REMOVE_DUPLICATES=false \
			TMP_DIR=${tmp_dirname}
		"""
}

process picard_index {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_SAMTOOLS

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir PICARD_DUPLICATE_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_bam_index

	output:
		set val(sample), file(bam), file("*.bai") \
			into \
				picard_bai_transcript_integrity,
				picard_bai_rnaseqc,
				picard_bai_tss,
				picard_bai_stats,
				picard_bai_idxstats,
				picard_bai_multiqc

	shell:
		"""
		samtools index ${bam}
		"""
}

process picard_complexity {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_PICARD

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir PICARD_COMPLEXITY_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_bam_complexity

	output:
		set val(sample), file("*.complexity") into picard_complexity

	shell:

		// the temporary dictory and the metrics file
		tmp_dirname = "tmp"
		metrics_filename = sample + ".complexity"

		"""
		java -Xmx10g -Djava.io.tmpdir=${tmp_dirname} \
			-jar \$EBROOTPICARD/picard.jar EstimateLibraryComplexity \
			VALIDATION_STRINGENCY=SILENT \
			INPUT=${bam} \
			OUTPUT=${metrics_filename} \
			TMP_DIR=${tmp_dirname}
		"""
}

process picard_rnaseqmetrics {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_PICARD

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir PICARD_RNASEQMETRICS_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_bam_rnaseqmetrics

	output:
		set val(sample), file("*.rnaseqmetrics") into picard_rnaseqmetrics

	shell:

		// the temporary dictory and the metrics file
		tmp_dirname = "tmp"
		metrics_filename = sample + ".rnaseqmetrics"

		"""
		java -Xmx10g -Djava.io.tmpdir=${tmp_dirname} \
			-jar \$EBROOTPICARD/picard.jar CollectRnaSeqMetrics \
			VALIDATION_STRINGENCY=SILENT \
			INPUT=${bam} \
			OUTPUT=${metrics_filename} \
			REF_FLAT=${ANNOT_REFFLAT_FILEPATH} \
			STRAND_SPECIFICITY=NONE \
			RIBOSOMAL_INTERVALS=${ANNOT_RRNA_INTERVAL_FILEPATH} \
			TMP_DIR=${tmp_dirname}
		"""
}

process picard_multimetrics {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_PICARD
	module MODULE_R

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir PICARD_MULTIMETRICS_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_bam_multimetrics

	output:
		set val(sample), file("*.pdf") into picard_multimetrics_pdf
		set val(sample), file("*_metrics") into picard_multimetrics_metrics

	shell:

		// the temporary dictory and the metrics file
		tmp_dirname = "tmp"
		metrics_filename = sample + ".multimetrics"

		"""
		java -Xmx10g -Djava.io.tmpdir=${tmp_dirname} \
			-jar \$EBROOTPICARD/picard.jar CollectMultipleMetrics \
			VALIDATION_STRINGENCY=SILENT \
			INPUT=${bam} \
			OUTPUT=${metrics_filename} \
			PROGRAM=CollectAlignmentSummaryMetrics \
			R=${SEQ_FILEPATH} \
			TMP_DIR=${tmp_dirname}
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             SAMTOOLS                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process stats {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_SAMTOOLS

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir STATS_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam), file(bai) from picard_bai_stats

	output:
		set val(sample), file("*.stats") into stats

	shell:
		"""
		samtools stats ${bam} > ${sample}.stats
		"""
}

process idxstats {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_SAMTOOLS

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir IDXSTATS_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam), file(bai) from picard_bai_idxstats

	output:
		set val(sample), file("*.idxstats") \
			into \
				idxstats_chrom,
				idxstats_multiqc

	shell:
		"""
		samtools idxstats ${bam} > ${sample}.idxstats
		"""
}

process chromosome {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_ANACONDA
	conda ANACONDA_ENV

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir CHROM_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(idxstats) from idxstats_chrom

	output:
		set val(sample), file("*_chrom.json") into chrom
	
	shell:
		"""
		python ${CHROM_SCRIPT} ${idxstats} > ${sample}_chrom.json
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                              RSEQC                                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process infer_experiment {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir RSEQC_INFER_EXPERIMENT_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_bam_infer_experiment
	
	output:
		set val(sample), file("*.infer_experiment*") into infer_experiment
	
	shell:

		metrics_filename = sample + ".infer_experiment"

		"""
		infer_experiment.py \
			-i ${bam} \
			-r ${ANNOT_BED_FILEPATH} \
			> ${metrics_filename}
		"""
}

process junction_annotation {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir RSEQC_JUNCTION_ANNOTATION_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_bam_junction_annotation

	output:
		set val(sample), file("*.junction_annotation*") into junction_annotation
	
	shell:

		metrics_filename = sample + ".junction_annotation"

		"""
		junction_annotation.py \
			-i ${bam} \
			-r ${ANNOT_BED_FILEPATH} \
			-o ${metrics_filename}
		"""
}

process junction_saturation {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir RSEQC_JUNCTION_SATURATION_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_bam_junction_saturation

	output:
		set val(sample), file("*.junction_saturation*") into junction_saturation

	shell:

		metrics_filename = sample + ".junction_saturation"

		"""
		junction_saturation.py \
			-i ${bam} \
			-r ${ANNOT_BED_FILEPATH} \
			-o ${metrics_filename}
		"""
}

process mismatch_profile {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir RSEQC_MISMATCH_PROFILE_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_bam_mismatch_profile

	output:
		set val(sample), file("*.mismatch_profile*") into mismatch_profile

	shell:

		metrics_filename = sample + ".mismatch_profile"

		"""
		mismatch_profile.py \
			-i ${bam} \
			-l ${ROUGH_READ_LENGTH} \
			-o ${metrics_filename}
		"""
}

process read_distribution {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir RSEQC_READ_DISTRIBUTION_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_bam_read_distribution

	output:
		set val(sample), file("*.read_distribution*") into read_distribution
	
	shell:

		metrics_filename = sample + ".read_distribution"

		"""
		read_distribution.py \
			-i ${bam} \
			-r ${ANNOT_BED_FILEPATH} \
			> ${metrics_filename}
		"""
}

process transcript_integrity {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir RSEQC_TRANSCRIPT_INTEGRITY_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam),
			file(bai) from picard_bai_transcript_integrity

	output:
		set val(sample), file("*.{xls,txt}") into transcript_integrity

	shell:
		"""
		tin.py \
			-i ${bam} \
			-r ${ANNOT_BED_FILEPATH}
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             RNASEQC                                 -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process rnaseqc {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RNASEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir RNASEQC_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam), file(bai) from picard_bai_rnaseqc

	output:
		set val(sample), file("*rnaseqc*") into rnaseqc

	shell:
		
		metrics_filename = sample + ".rnaseqc"

		if (SINGLE_END) {

			single_end = "-singleEnd"

		} else {

			single_end = ""

		}

		"""
		java -Xmx10g -jar \$EBROOTRNAMINSEQC/RNA-SeQC_v[0-9].[0-9].[0-9].jar \
			-d 1000000 \
			-rRNA ${ANNOT_RRNA_FILEPATH} \
			-r ${SEQ_FILEPATH} \
			-t ${ANNOT_RNASEQC_FILEPATH} \
			-o ${metrics_filename} \
			-gatkFlags '-S SILENT -U ALLOW_SEQ_DICT_INCOMPATIBILITY' \
			${single_end} \
			-s "${sample}|${bam}|${sample}"
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                           TSS COVERAGE                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process tss {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_ANACONDA
	conda ANACONDA_ENV

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir TSS_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE
	
	input:
		set val(sample), file(bam), file(bai) from picard_bai_tss
	
	output:
		set val(sample), file("*_tss.json") into tss
	
	shell:
		"""
		python ${TSS_SCRIPT} \
			--bam ${bam} \
			--gtf ${ANNOT_GTF_FILEPATH} \
			--win 3000 \
			--size 200 \
			--name ${sample} \
			> ${sample}_tss.json
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       QC GENES STRANDEDNESS                         -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process qc_genes {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_GCC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir QC_GENES_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE
	
	input:
		file info from info_qc_genes
		set val(sample), file(bam) from picard_bam_qc_genes
	
	output:
		set val(sample), file("*_qc_genes.json") into qc_genes
	
	shell:
		"""
		${QC_GENES_SCRIPT} \
			${bam} \
			${ANNOT_BED_FILEPATH} \
			${RRNA} \
			${MTRRNA} \
			${GLOBIN} \
			> ${sample}_qc_genes.json
		"""
}

process transcripts {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_GCC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir TRANSCRIPTS_DIRPATH,
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE
	
	input:
		file info from info_transcripts
		set val(sample), file(bam) from picard_bam_transcripts
	
	output:
		set val(sample), file("*_transcripts.json") into transcripts
	
	shell:
		"""
		${TRANSCRIPTS_SCRIPT} \
			${bam} \
			${ANNOT_BED_FILEPATH} \
			> ${sample}_transcripts.json
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                            MULTIQC                                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process multiqc_conf {

	// problem with slurm otherwise
	beforeScript "module purge"

	// module
	module MODULE_PYTHON

	// HPC
	cpus 1
	executor "slurm"
	memory "1000"

	output:
		file "multiqc_conf.yml" into multiqc_conf

	shell:

		def filename = "multiqc_conf.yml"

		def logo = Paths.get(
						workflow.projectDir.toString(),
						"conf",
						"logo.jpg"
						).toString()

		def h = [:]
		h["Contact e-mail"] = "bioinformatics@crick.ac.uk"
		h["Experiment type"] = "RNA-Seq"
		h["Protocol"] = SINGLE_END ? "single-end" : "paired-end"
		h["Parametrised species"] = SPECIES
		h["Parametrised strandedness"] = STRANDEDNESS
		h["Parametrised directory"] = DIRECTORY
		h["Determined read length"] = READ_LENGTH
		//h["Inferred strandedness"] = INFERRED_STRANDEDNESS
		h["GTF general annotation file"] = ANNOT_GTF_FILEPATH
		h["FASTA genome sequence file"] = SEQ_FILEPATH
		h["RSEM-STAR index files prefix"] = RSEM_STAR_INDICE_PREFIX_FILEPATH
		h["BED genome intervals file"] = ANNOT_BED_FILEPATH
		h["REFFLAT gene annotation file"] = ANNOT_REFFLAT_FILEPATH
		h["GTF RNASeQC specific file"] = ANNOT_RNASEQC_FILEPATH
		h["Ribosomal RNAs loci file"] = ANNOT_RRNA_FILEPATH
		h["Ribosomal RNAs intervals file"] = ANNOT_RRNA_INTERVAL_FILEPATH


		def conf = [:]
		conf["custom_logo"] = logo
		conf["report_header_info"] = h

		def json = JsonOutput.toJson(conf)

		// FAIL TO CONVERT NESTED MAP/DICT
		//"""
		//echo -E '${json}' > multiqc_conf.json
		//python -c '
		//import sys, yaml, json
		//yaml.safe_dump(
		//	json.load(sys.stdin),
		//	sys.stdout,
		//	default_flow_style=False
		//	)
		//' < multiqc_conf.json > conf.yml
		//cp -v ${MULTIQC_CONF_FILEPATH} ${filename}
		//cat conf.yml >> ${filename}
		//"""

		def header = []
		h.each{ k, v ->
			s = "   - " + k + ": " + v
			header.add(s)
		}
		def header_string = header.join("\\n")

		"""
		cp -v ${MULTIQC_CONF_FILEPATH} ${filename}
		echo "custom_logo: ${logo}" >> ${filename}
		echo "report_header_info:" >> ${filename}
		echo -e "${header_string}" >> ${filename}
		"""
}

process multiqc {

	echo true

	// problem with slurm otherwise
	beforeScript "module purge"

	// module
	//module MODULE_PYTHON
	//module MODULE_MULTIQC
	module MODULE_ANACONDA

	// anaconda
	conda ANACONDA_ENV

	// HPC
	cpus 1
	executor "slurm"
	memory "6000"

	// output directory
	publishDir workflow.projectDir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		
		file (conf) from multiqc_conf
		
		file("*") from cutadapt_log.collect()
		file("*") from cutadapt_multiqc.map{x->x[1]}.collect()

		file("*") from fastqc.collect()
		file("*") from fastqc_cutadapt.collect()

		file("*") from strandedness_multiqc.collect()

		file("*") from fscreen_html.map{x->x[1]}.collect()
		file("*") from fscreen_txt.map{x->x[1]}.collect()

		file("*") from rsem_transcript.map{x->x[1]}.collect()
		file("*") from rsem_results.map{x->x[1]}.collect()
		file("*") from rsem_stat.map{x->x[1]}.collect()
		file("*") from rsem_star_genome_multiqc.map{x->x[1]}.collect()

		file("*") from picard_rg_multiqc.map{x->x[1]}.collect()
		file("*") from picard_duplicate.map{x->x[1]}.collect()
		file("*") from picard_bai_multiqc.map{x->x[1]}.collect()
		file("*") from picard_complexity.map{x->x[1]}.collect()
		file("*") from picard_rnaseqmetrics.map{x->x[1]}.collect()
		file("*") from picard_multimetrics_pdf.map{x->x[1]}.collect()
		file("*") from picard_multimetrics_metrics.map{x->x[1]}.collect()

		file("*") from stats.collect()
		file("*") from idxstats_multiqc.collect()

		file("*") from infer_experiment.map{x->x[1]}.collect()
		file("*") from junction_annotation.map{x->x[1]}.collect()
		file("*") from junction_saturation.map{x->x[1]}.collect()
		file("*") from mismatch_profile.map{x->x[1]}.collect()
		file("*") from read_distribution.map{x->x[1]}.collect()
		file("*") from transcript_integrity.map{x->x[1]}.collect()

		file("*") from rnaseqc.map{x->x[1]}.collect()

		file("*") from chrom.map{x->x[1]}.collect()
		file("*") from tss.map{x->x[1]}.collect()
		file("*") from info_qc_genes_multiqc.collect()
		file("*") from qc_genes.map{x->x[1]}.collect()
		file("*") from info_transcripts_multiqc.collect()
		file("*") from transcripts.map{x->x[1]}.collect()

	output:
		file "multiqc_data" into multiqc_data
		file "multiqc_report.html" into multiqc_report

	shell:
		"""
		export LC_ALL=en_US.utf8
		export LANG=en_US.utf8
		multiqc \
			--data-format json \
			--template custom \
			--config ${conf} \
			.
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                           NOTIFICATION                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

workflow.onComplete {

	// the bioinformatician
	def USER = System.getenv("USER")

	// might need to be changed but don't want to give name of the
	// bioinformatician as a parameter
	def recipient = USER + "@crick.ac.uk"

	// if it failed or not
	def subject = "[nexflow] the QC workflow has successfuly finished"
	if (!workflow.success) {
		subject = "[nextflow] PROBLEM: the QC workflow has failed"
	}

	// just the strict necessary, need to add more
	def body =
	"""
	Start time: ${workflow.start}
	End time: ${workflow.complete}
	Exit status: ${workflow.exitStatus}
	"""

	// if it failed... when it fails
	if (workflow.errorMessage) {
		body = body +
	"""
	Error message: ${workflow.errorMessage}
	"""
	}

	// if it failed... when it fails
	if (workflow.errorReport) {
		body = body +
	"""
	Report message: ${workflow.errorReport}
	"""
	}
	
	// there is only mutt installed on the node
	["mutt", "-s", subject, recipient].execute() << body
}

