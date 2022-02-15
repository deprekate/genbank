# genbank
Python code to work with Genbank files

This repo contains several classes to help work with Genbank files

The flow goes:
```
File > Read > Locus > Feature
```

To use:
```
from genbank.file import File

f = File('infile.gbk')
for name,locus in f.items():
	print("-----Locus:",name, "-----")
	for feature in locus:
		print(feature)
```
