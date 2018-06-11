#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
		name = "multiqc_custom",
		version = "0.1",
		author = "nourdine bah",
		author_email = "nourdine.bah@crick.ac.uk",
		description = "MultiQC plugin for special metrics",
		packages = find_packages(),
		include_package_data = True,
		install_requires = ["multiqc==1.5"],
		entry_points = {
			"multiqc.modules.v1": [
				"tss = qc_plugin.modules.tss:MultiqcModule",
				"qc_genes = qc_plugin.modules.qc_genes:MultiqcModule",
				"chrom = qc_plugin.modules.chrom:MultiqcModule",
				"transcripts = qc_plugin.modules.transcripts:MultiqcModule"
				],
			"multiqc.templates.v1": [
				"custom = qc_plugin.templates.custom"
				]
			}
		)
