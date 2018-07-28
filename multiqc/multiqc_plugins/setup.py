#!/usr/bin/env python

from setuptools import setup, find_packages

html_files = ["multiqc_plugins/templates/custom/includes.html"]
js_files = \
		["multiqc_plugins/templates/custom/assets/js/highcharts-more.v5.0.6.js"]

setup(
		name = "multiqc_plugins",
		version = "0.1",
		author = "nourdine bah",
		author_email = "nourdine.bah@crick.ac.uk",
		description = "MultiQC plugins for special metrics",
		packages = find_packages(),
		include_package_data = True,
		install_requires = ["multiqc>=1.5"],
		data_files=[("html", html_files), ("js", js_files)],
		entry_points = {
			"multiqc.modules.v1": [
				"tss = multiqc_plugins.modules.tss:MultiqcModule",
				"qc_genes = multiqc_plugins.modules.qc_genes:MultiqcModule",
				"chrom = multiqc_plugins.modules.chrom:MultiqcModule",
				"transcripts = multiqc_plugins.modules.transcripts:MultiqcModule"
				],
			"multiqc.templates.v1": [
				"custom = multiqc_plugins.templates.custom"
				]
			}
		)

