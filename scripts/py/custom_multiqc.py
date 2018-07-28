from multiqc.plots import bargraph

data = {
		'sample 1': {
			'aligned': 23542,
			'not_aligned': 343,
			},
		'sample 2': {
			'not_aligned': 7328,
			'aligned': 1275,
			}
		}
html_content = bargraph.plot(data)
print(html_content)
