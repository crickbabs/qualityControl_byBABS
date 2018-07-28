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

### If you are from BABS

You don't need to do anything. The pipeline is using a `qc_pipeline` conda  environment available in our shared space : `/camp/stp/babs/working/software/anaconda/envs`. Just be sure that anaconda is properly configured.

	$ readlink $HOME/.conda
	/camp/stp/babs/working/$USER/.conda

	$ readlink $HOME/.condarc
	/camp/stp/babs/working/$USER/.condarc

	$ cat /camp/stp/babs/working/$USER/.condarc
	envs_dirs:
	 - /camp/stp/babs/working/software/anaconda/envs
	pkgs_dirs:
	 - /camp/stp/babs/working/software/anaconda/pkgs

	$ cat /camp/stp/babs/working/$USER/.conda/environments.txt
	...
	/camp/stp/babs/working/software/anaconda/envs/qc_pipeline
	...

### If you are not from BABS

You need to create an anaconda environment and install multiqc, other packages, and the multiqc plugins inside this environment. Here is the procedure.

First, you need to load Anaconda:

	$ module load Anaconda2/5.1.0

Then, you create a new environment named `qc_pipeline`:

	$ conda create --yes --name qc_pipeline python=3.6

After that you can install somes packages that the pipeline requires:

	$ conda install --yes --channel bioconda --name qc_pipeline multiqc=1.5
	$ conda install --yes --channel anaconda --name qc_pipeline openblas=0.2.20
	$ conda install --yes --channel bioconda --name qc_pipeline htseq=0.9.1
	$ conda install --yes --channel anaconda --name qc_pipeline pandas=0.23.3

Install the multiqc plugins:

	$ conda install --yes --name qc_pipeline multiqc/multiqc_plugins-1.0-py36_2.tar.bz2

Finally the variable `ANACONDA_ENV` needs to be changed in the `main.nf` file. This variable has to be set to the path of the anaconda environment you just created.

## Configuring the pipeline

For example:

	$ echo "directory: $(realpath fastq/paired_end/)" >> params.yml

## Running the pipeline
	
	$ module load nextflow/0.30.0
	$ nextflow run -params-file params.yml main.nf -with-timeline

