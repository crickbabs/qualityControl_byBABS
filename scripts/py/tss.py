#!/usr/bin/python

import sys
from os import listdir
from os.path import join
import json
import pandas as pd


JSON_DIRPATH = "/camp/stp/babs/working/bahn/code/nf/qc_pipeline/output/tss"

# get the json files
files = sorted([ join(JSON_DIRPATH, f) for f in listdir(JSON_DIRPATH) ])
jsons = list(map( lambda x: json.load(open(x)), files ))

# merge all the samples in one dict
data = {}
for j in jsons:
	d = {}
	for x, y in zip( j["support"], j["coverage"] ):
		d[int(x)] = int(y)
	data[ j["name"] ] = d

# output as a json
output = {
		"id": "tss_lineplot",
		"section_name": "TSS coverage",
		"description": "The coverage of the TSS region for each sample.",
		"plot_type": "linegraph",
		"pconfig": {
			"id": "tss_linegraph",
			"title": "TSS coverage",
			"ylab": "Coverate (bp)"
			},
		"data": data
		}
print( json.dumps(output, indent=3) )

