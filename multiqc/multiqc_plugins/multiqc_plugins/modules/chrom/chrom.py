#!/usr/bin/env python

# nourdine.bah@crick.ac.uk

from __future__ import print_function
from collections import OrderedDict
import logging
import json
#from unicodedata import normalize

from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# the multiqc logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):

	def __init__(self):

		super(MultiqcModule, self).__init__(
				name = "Chromosome representativity",
				target = "",
				anchor = "chrom",
				href = "",
				info = ""
				)

		self.data = list()
		
		#########################################################################
		for f in self.find_log_files("chrom"):
			d = self.parse_json(f)
			self.data.append(d)
			######################################################################

		if len(self.data) == 0:
			log.debug(
					"Could not find any reports in {}".format(config.analysis_dir)
					)
			raise UserWarning

		log.info("Found {} reports".format(len(self.data)))
		#########################################################################

		# create highcharts plot categories and series
		self.categories = [ "chr" + c for c in self.data[0]["sequence_name"] ]
		self.series = list()
		for s in self.data:
			d = dict()
			d["name"] = s["name"]
			d["data"] = s["diff"]
			self.series.append(d)

		#########################################################################
		# PLOT

		description = """For each chromosome, the difference between the
		proportion of mapped reads and the expected value of this metrics based
		on the chromosome length with an uniform distribution."""

		mailto = """<a href="mailto:bioinformatics@crick.ac.uk">
		bioinformatics@crick.ac.uk</a>"""

		self.add_section(
				description = description,
				helptext = mailto,
				plot = self.chrom_plot()
				)

	############################################################################
	def parse_json(self, f):
		j = json.loads(f['f'])
		j["name"] = f["fn"].replace("_chrom.json", "")
		return j

	############################################################################
	def chrom_plot(self):
		
		html = """
		<div id="chrom_plot" class="hc-plot"></div>
		<script type="text/javascript">
			var chrom_categories = """ + json.dumps(self.categories) + """;
			var chrom_series = """ + json.dumps(self.series) + """;
			$(function () {
				$("#chrom_plot").highcharts({
					title: { text:
					"Deviation of the reads proportion from its expected value"
					},
					chart: { type: "column" },
					yAxis: {
						title: { text: "Deviation of the reads proportion" }
					},
					xAxis: {
						title: { text: "Chromosome" },
						categories: chrom_categories
					},
					series: chrom_series
				});
			});
		</script>
		"""

		return html

	############################################################################
	def test_plot(self):
		
		html = """
		<div id="chrom_test_plot" class="hc-plot"></div>
		<script type="text/javascript">
			var chrom_categories = ["chrom1", "chrom2", "chrom3", "chrom4"];
			var chrom_series = [
				{name:"sample1", data:[1, 3, -6, 2]},
				{name:"sample2", data:[4, 6, -1, 6]},
				{name:"sample3", data:[4, 3, -9, 9]},
				{name:"sample4", data:[4, 2, -3, 6]}
			];
			$(function () {
				$("#chrom_test_plot").highcharts({
					chart: { type: "column" },
					xAxis: { categories: chrom_categories },
					series: chrom_series
				});
			});
		</script>
		"""

		return html

