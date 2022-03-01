import sys

from genbank.file import File

f = File(sys.argv[1])
print(f)
for name,locus in f.items():
	print("-----Locus:", name, "-----")
	for feature in locus:
		print(feature)

	print(locus.nearest(59, -1, 'aaa'))
	print(locus.distance(59, -1, 'aaa'))
