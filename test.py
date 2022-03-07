import sys

from genbank.file import File

f = File(sys.argv[1])
for name,locus in f.items():
	print("-----Locus:", name, "-----")
	for feature in locus:
		print(feature)
		'''
		if feature.type_is('CDS'):
			for codon in feature.codons():
				print(codon)
		'''
	stops = ['taa','tag','tga']
	print(locus.last(59, +1, stops))
	print(locus.last(59, +1, 'taa'))

