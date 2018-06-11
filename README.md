# qualityControl-byBABS

## Obtaining the Files You Need

	$ git clone https://github.com/crickbabs/qualityControl_byBABS
	$ cd qualityControl_byBABS

## Compiling two required binaries
	
	$ cd script/cpp
	$ sh qc_directionality.sh compil
	$ sh transcript_directionality.sh compil
	$ cd -

## Installing the custom MultiQC module

	$ cd script/multiqc/qc-plugin
	$ module load Python/3.6.3-foss-2017b
	$ python setup.py develop --user
	$ cd -

## Configuring the pipeline

For example:

	$ echo "directory: $(realpath fastq/paired_end/)" >> params.yml

## Running the pipeline
	
	$ module load nextflow/0.30.0
	$ nextflow run -params-file params.yml -timeline main.nf

