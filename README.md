# qualityControl-byBABS

## Obtaining the Files You Need

	$ git clone https://github.com/crickbabs/qualityControl_byBABS
	$ cd qualityControl_byBABS

## Compiling two required binaries
	
	$ cd script/cpp
	$ sh qc_directionality.sh
	$ sh transcript_directionality.sh
	$ cd -

## Installing the custom MultiQC module

	$ cd script/multiqc/qc-plugin
	$ module load Python/3.6.3-foss-2017b
	$ python setup.py develop --user
	$ cd -

## Configuring the pipeline

	$ vim params.yml

## Running the pipeline
	
	$ module load nextflow/0.30.0
	$ nextflow run -params-file params.yml main.nf

