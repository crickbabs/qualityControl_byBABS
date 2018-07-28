#!/usr/bin/env python

# nourdine.bah@crick.ac.uk

from __future__ import print_function
from collections import OrderedDict
import logging
import json
import numpy as np

import random
from random import randrange

from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# the multiqc logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):

	def __init__(self):

		super(MultiqcModule, self).__init__(
				name = "Reads strandedness per transcript",
				target = "",
				anchor = "transcripts",
				href = "",
				info = ""
				)

		self.data = list()
		
		#########################################################################
		for f in self.find_log_files("transcripts"):

			if f["fn"] == "experiment_info_transcripts.json":
				info = json.loads(f['f'])
				continue

			d = self.parse_json(f)

			######################################################################
			## test
			#dd = {}
			#for k, v in d["data"].items():
			#	x = round(random.gauss(mu=10, sigma=10), 3)
			#	ddd = {"sense": x, "antisense": 1.0-x, "undetermined": 0.0}
			#	dd[k] = ddd
			#d["data"] = dd
			######################################################################

			self.data.append(d)
			######################################################################

		if len(self.data) == 0:
			log.debug(
					"Could not find any reports in {}".format(config.analysis_dir)
					)
			raise UserWarning

		log.info("Found {} reports".format(len(self.data)))
		#########################################################################

		## SCIENTIST DEFINED STRANDEDNESS (not the one guessed from the pipeline)

		if info["strandedness"].lower() == "forward":
			self.good_strandedness = "sense"
			self.bad_strandedness = "antisense"
			self.protocol = "forward"
			self.direction = self.good_strandedness
			self.value = "1"
		elif info["strandedness"].lower() == "reverse":
			self.good_strandedness = "antisense"
			self.bad_strandedness = "sense"
			self.protocol = "reverse"
			self.direction = self.good_strandedness
			self.value = "1"
		else:
			self.good_strandedness = "sense"
			self.bad_strandedness = "antisense"
			self.protocol = "none"
			self.direction = self.good_strandedness
			self.value = "0.5"

		#########################################################################

		self.plot_data = list()

		for sample in self.data:

			n = sample["name"]
			d = sample["data"]

			# storage
			dd = dict()
			data = dict()
			inter_quartile_range = dict()
			upper_thresh = dict()
			lower_thresh = dict()
			boxplot = dict()
			outliers = dict()

			# merge the values
			data["sense"] = [ v["sense"] for k, v in d.items() ]
			data["antisense"] = [ v["antisense"] for k, v in d.items() ]
			data["undetermined"] = [ v["undetermined"] for k, v in d.items() ]

			# threshold for the outliers
			for k, v in data.items():

				boxplot[k] = list()
				outliers[k] = list()

				iqr = np.percentile(v, 75) - np.percentile(v, 25)
				lt = np.percentile(v, 25) - iqr * 1.5
				ut = np.percentile(v, 75) + iqr * 1.5

				for x in v:
					if x < lt or ut < x:
						outliers[k].append(x)
					else:
						boxplot[k].append(x)

				inter_quartile_range[k] = iqr
				lower_thresh[k] = lt
				upper_thresh[k] = ut

			dd["name"] = n
			dd["boxplot"] = boxplot[self.good_strandedness]

			#dd["outliers"] = outliers[self.good_strandedness]
			dd["outliers"] = []
			m = 100
			if len(outliers[self.good_strandedness]) > m:
				for i in range(m):
					n = randrange( len(outliers[self.good_strandedness]) )
					dd["outliers"].append(outliers[self.good_strandedness][n])

			self.plot_data.append(dd)

		#########################################################################

		self.plot_data = sorted(self.plot_data, key=lambda x: x["name"])
		self.categories = [ x["name"] for x in self.plot_data ]
		self.boxplot = list()
		# !!! NEED TO ADD THE X INDEX AT THE BEGINING OF THE ARRAY
		for i, x in enumerate(self.plot_data):
			self.boxplot.append( [i] + sorted(x["boxplot"]) )
		self.outliers = list()
		for i, d in enumerate(self.plot_data):
			for x in d["outliers"]:
				self.outliers.append( [i, x] )

		#########################################################################
		# PLOT

		description = """Boxplots of the proportion of the <b>{d}</b> reads
		per mapped transcript. The protocol of the experiment is supposed to be
		<b>{p}</b>, so we can expect to have the values around <b>{v}</b>. Only
		<b>1000 outliers per sample</b> have been randomly selected and plotted
		here.
		""".format(**{"d": self.direction, "v": self.value, "p": self.protocol})

		mailto = """<a href="mailto:bioinformatics@crick.ac.uk">
		bioinformatics@crick.ac.uk</a>"""

		self.add_section(
				description = description,
				helptext = mailto,
				plot = self.transcripts_plot()
				)

	############################################################################
	def parse_json(self, f):
		d = dict()
		d["name"] = f["fn"].replace("_transcripts.json", "")
		d["data"] = json.loads(f['f'])
		return d

	############################################################################
	def transcripts_plot(self):
		
		html = """
		<div id="transcripts_plot" class="hc-plot"></div>
		<script type="text/javascript">
			var transcripts_categories = """ + json.dumps(self.categories) + """;
			var transcripts_data = """ + json.dumps(self.boxplot) + """;
			var transcripts_outliers = """ + json.dumps(self.outliers) + """;
			$(function () {
				$("#transcripts_plot").highcharts({
					title: { text: "Reads strandedness" },
					yAxis: {
						title: {
							text: "Proportion of """ + self.direction + """ reads"
							}
					},
					xAxis: {
						title: { text: "Sample" },
						categories: transcripts_categories
					},
					series: [
					{
						name: "Boxplot",
						type: "boxplot",
						data: transcripts_data
					},
					{
						name: "Outlier",
						type: "scatter",
						data: transcripts_outliers
					}]
				});
			});
		</script>
		"""

		return html

	############################################################################
	def test_plot(self):
		
		html = """
		<div id="transcripts_test_plot" class="hc-plot"></div>
		<script type="text/javascript">
			$(function () {
				$("#transcripts_test_plot").highcharts({
					xAxis: {
						categories: ["sample1", "sample2", "sample3", "sample4"]
						},
					series: [
					{
						name: "Observations",
						type: "boxplot",
						data: [
						[0, 760, 801, 848, 895, 965],
						[1, 733, 853, 939, 980, 1080],
						[2, 724, 802, 806, 871, 950],
						[3, 834, 836, 864, 882, 910]
						]
					},
					{
						name: "Outlier",
						type: "scatter",
						data: [ [0, 644], [3, 718], [3, 951], [3, 969] ]
					}]
				});
			});
		</script>
		"""

		return html

