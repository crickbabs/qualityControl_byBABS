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
				name = "Strandedness of QC genes",
				target = "",
				anchor = "qc_genes",
				href = "",
				info = ""
				)

		self.data = list()
		
		#########################################################################
		for f in self.find_log_files("qc_genes"):

			if f["fn"] == "experiment_info_qc_genes.json":
				info = json.loads(f['f'])
				continue

			d = self.parse_json(f)

			# just compute the number of unidentified reads
			d["unidentified"] = d["mapped"] - d["success"]

			# compute the total strandedness
			d["sense"] = \
					d["sense globin"] +\
					d["sense rRNA"] +\
					d["sense mt-rRNA"] +\
					d["sense non QC"]
			d["antisense"] = \
					d["antisense globin"] +\
					d["antisense rRNA"] +\
					d["antisense mt-rRNA"] +\
					d["antisense non QC"]

			self.data.append(d)
			######################################################################

		if len(self.data) == 0:
			log.debug(
					"Could not find any reports in {}".format(config.analysis_dir)
					)
			raise UserWarning

		log.info("Found {} reports".format(len(self.data)))
		#########################################################################

		## QC GENES STRANDEDNESS

		# sort by the column size
		self.data = sorted(self.data, key= lambda x: x["mapped"])

		# reformat the data for plotting, to understand that you can check
		# the self.test_plot method
		self.categories = list(map( lambda x: x["name"], self.data ))
		self.series = list()
		metrics = [
				"sense globin",
				"antisense globin",
				"sense mt-rRNA",
				"antisense mt-rRNA",
				"sense rRNA",
				"antisense rRNA",
				"sense non QC",
				"antisense non QC",
				"unidentified",
				"total"
				]
		metrics.reverse()
		for m in metrics:
			plot_type = "scatter" if m == "total" else "column"
			data = list(map( lambda x: x[m], self.data ))
			d = {"type": plot_type, "name": m, "data": data}
			self.series.append(d)

		#########################################################################
		# PLOT

		description = """The strandedness of the reads mapping to some quality
		control genes (ribosomal RNAs, mitochondrial ribosomal RNAs and globin).
		"""

		mailto = """<a href="mailto:bioinformatics@crick.ac.uk">
		bioinformatics@crick.ac.uk</a>"""

		self.add_section(
				description = description,
				helptext = mailto,
				plot = self.qc_genes_plot()
				)

	############################################################################
	def parse_json(self, f):
		j = json.loads(f['f'])
		j["name"] = f["fn"].replace("_qc_genes.json", "")
		return j

	############################################################################
	def qc_genes_plot(self):
		
		html = """
		<div id="qc_genes_plot" class="hc-plot"></div>
		<script type="text/javascript">
			var qc_genes_categories = """ + json.dumps(self.categories) + """;
			var qc_genes_series = """ + json.dumps(self.series) + """;
			$(function () {
				$("#qc_genes_plot").highcharts({
					title: { text: "Strandedness" },
					yAxis: { title: { text: "Read count" } },
					xAxis: {
						title: { text: "Sample" },
						categories: qc_genes_categories
					},
					plotOptions: {
						column: { stacking: "normal" }
					},
					series: qc_genes_series
				});
			});
		</script>
		"""

		return html

	############################################################################
	def qc_genes_test_plot(self):
		
		html = """
		<div id="qc_genes_test_plot" class="hc-plot"></div>
		<script type="text/javascript">
			$(function () {
				$("#qc_genes_test_plot").highcharts({
					xAxis: {
						categories: ["sample1", "sample2", "sample3", "sample4"]
						},
					plotOptions: {
						column: { stacking: "normal" }
					},
					series: [
					{ type: "column", name: "sense rRNA", data: [5,3,2,4] },
					{ type: "column", name: "antisense rRNA", data: [8,4,2,2] },
					{ type: "column", name: "sense globin", data: [9,4,3,8] },
					{ type: "column", name: "antisense globin", data: [7,5,6,9] },
					{ type: "scatter", "name": "total", data: [20, 29, 24, 26] },
					{ type: "scatter", "name": "success", data: [15, 16, 13, 19] }
					]
				});
			});
		</script>
		"""

		return html

	############################################################################
	def total_test_plot(self):
		
		html = """
		<div id="total_test_plot" class="hc-plot"></div>
		<script type="text/javascript">
			$(function () {
				$("#total_test_plot").highcharts({
					xAxis: {
						categories: ["sample1", "sample2", "sample3", "sample4"]
						},
					yAxis: { max: 1 },
					plotOptions: {
						column: { stacking: "normal" }
					},
					series: [
					{ type: "column", name: "sense", data: [.5,.4,.1] },
					{ type: "column", name: "antisense", data: [.2,.2,.6] },
					{ type: "column", name: "undetermined", data: [.3,.4,.3] }
					]
				});
			});
		</script>
		"""

		return html

