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




This package also allows you to perform various operations on a given genome file:
```
$ genbank.py tests/phiX174.gbk -f tabular
'phiX174'	'CDS'	(('100', '627'),)	{'gene': '"G"'}	a
'phiX174'	'CDS'	(('636', '1622'),)	{'gene': '"H"'}
'phiX174'	'CDS'	(('1659', '3227'),)	{'gene': '"A"'}
'phiX174'	'CDS'	(('2780', '3142'),)	{'gene': '"B"'}
'phiX174'	'CDS'	(('3142', '3312'),)	{'gene': '"K"'}
...

$ genbank.py tests/phiX174.gbk -f fasta
>phiX174
gtgtgaggttataacgccgaagcggtaaaaattttaatttttgccgctgagggg
ttgaccaagcgaagcgcggtaggttttctgcttaggagtttaatcatgtttcag
...

$ genbank.py tests/phiX174.gbk -f fna
>phiX174_CDS_[100..627] [gene="G"]
atgtttcagacttttatttctcgccataattcaaactttttttctgataag
>phiX174_CDS_[636..1622] [gene="H"]
atgtttggtgctattgctggcggtattgcttctgctcttgctggtggcgcc
>phiX174_CDS_[1659..3227]
...

$ genbank.py tests/phiX174.gbk -f faa
>phiX174_CDS_[100..627] [gene="G"]
MFQTFISRHNSNFFSDKLVLTSVTPASSAPVLQTPKATSSTLYFDSLTVNA
>phiX174_CDS_[636..1622] [gene="H"]
MFGAIAGGIASALAGGAMSKLFGGGQKAASGGIQGDVLATDNNTVGMGDAG
>phiX174_CDS_[1659..3227] [gene="A"]
...
```
