#!/usr/bin/python

import json

rseqc_filepath = "/camp/stp/babs/working/bahn/projects/haydaya/maria.iannitto/test/paired_end_2/multiqc_data_1/multiqc_rseqc_infer_experiment.json"
rnaseqc_filepath = "/camp/stp/babs/working/bahn/projects/haydaya/maria.iannitto/test/paired_end_2/multiqc_data_1/multiqc_rna_seqc.json"

# load json as dict
with open(rseqc_filepath, 'r') as f:
	rseqc = json.load(f)
with open(rnaseqc_filepath, 'r') as f:
	rnaseqc = json.load(f)

# rseqc
rseqc_results = []
for k, v in rseqc.items():
	# we are getting the index with the alphabetical order (["failed",
	# "se/pe_antisense", "se/pe_sense"])
	keys = sorted(list( v.keys() ))
	sense = v[ keys[2] ]
	antisense = v[ keys[1] ]
	nrmsd_antisense = antisense / (sense + antisense) * 100.0
	nrmsd_sense = sense / (sense + antisense) * 100.0
	rseqc_results.append(nrmsd_sense)
rseqc_result = sum(rseqc_results) / len(rseqc_results)


# rnaseqc
rnaseqc_results = []
for k, v in rnaseqc.items():
	rnaseqc_results.append( (v["End 2 % Sense"] + v["End 2 % Sense"]) / 2.0 )
rnaseqc_result = sum(rnaseqc_results) / len(rnaseqc_results)






