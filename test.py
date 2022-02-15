import sys

from genbank.file import File

f = File(sys.argv[1])
for name,locus in f.items():
	print("-----Locus:",name, "-----")
	for feature in locus:
		print(feature)
