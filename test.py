import sys

from genbank.file import File

f = File(sys.argv[1])
print(f)
for name,locus in f.items():
	print("-----Locus:", name, "-----")
	for feature in locus:
		print(feature)

	stops = ('taa','tag','tga')
	print(locus.last(59, +1, stops))
	print(locus.last(59, +1, 'taa'))


