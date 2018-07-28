#!/usr/bin/env python

# nourdine.bah@crick.ac.uk

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# the multiqc logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):

	def __init__(self):

		super(MultiqcModule, self).__init__(
				name = "TSSs coverage",
				target = "",
				anchor = "tss",
				href = "",
				info = ""
				)

		self.tss_data = dict()
		self.coverage = list()
		self.support = list()
		
		#########################################################################
		for f in self.find_log_files("tss"):
			self.coverage.append( self.parse_coverage(f) )
			d = self.parse_json(f)
			self.support.append( d["support"] )
			######################################################################

		if len(self.coverage) == 0:
			log.debug(
					"Could not find any reports in {}".format(config.analysis_dir)
					)
			raise UserWarning

		log.info("Found {} reports".format(len(self.coverage)))

		#########################################################################
		# PLOT

		self.categories = self.support[0]

		description = """Just a plot of the transcription start sites coverage.
		The reads have been <b>lengthened by 200 base pairs</b>."""

		mailto = """<a href="mailto:bioinformatics@crick.ac.uk">
		bioinformatics@crick.ac.uk</a>"""

		self.add_section(
				description = description,
				helptext = mailto,
				plot = self.tss_plot()
				)

	############################################################################
	def parse_json(self, f):
		return json.loads(f['f'])

	############################################################################
	def parse_coverage(self, f):

		# f["fn"], f["s_name"]

		d = {}
		j = json.loads(f['f'])
		d["name"] = j["name"]
		d["data"] = j["coverage"]

		return d

	############################################################################
	def tss_plot(self):
		
		html = """
		<div id="tss_plot" class="hc-plot"></div>
		<script type="text/javascript">
			var tss_data = """ + json.dumps(self.coverage) + """;
			var tss_categories = """ + json.dumps(self.categories) + """;
			$(function () {
				$("#tss_plot").highcharts({
					type: "line",
					title: { text: "Transcription start sites coverage" },
					yAxis: { title: { text: "Base pair count" } },
					xAxis: {
						title: { text: "Base pair position" },
						categories: tss_categories
					},
					series: tss_data
				});
			});
		</script>
		"""

		return html

