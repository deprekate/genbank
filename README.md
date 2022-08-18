# genbank
Python code to work with Genbank files

This repo contains several classes to help work with Genbank files

The flow goes:
```
File -> Locus -> Feature
```

To use:
```python
from genbank.file import File

file = File('infile.gbk')
for locus in file:
	print(name)
	for feature in locus:
		print(feature)
```


You can also build a Locus object from the ground up:
```python
from genbank.locus import Locus
locus = Locus('test', 'actgactgatcgtagctagc')
# then add a feature by parsing text of a genbank feature
locus.read_feature('  CDS  1..10')
# or add one by manually specifing the type,strand,and locations
locus.add_feature('CDS',+1,[['10','20']])
locus.write()
```
which gives:
```
LOCUS       test        20 bp  DNA          UNK
DEFINITION  test
FEATURES          Location/Qualifiers
     CDS             1..10
     CDS             10..20
ORIGIN
        1 actgactgat cgtagctagc
//
```

---

This package also allows you to perform various conversions on a given genome file:
```bash
$ genbank.py tests/phiX174.gbk -f tabular
'phiX174'	'CDS'	(('100', '627'),)	{'gene': "G"}
'phiX174'	'CDS'	(('636', '1622'),)	{'gene': "H"}
'phiX174'	'CDS'	(('1659', '3227'),)	{'gene': "A"}
'phiX174'	'CDS'	(('2780', '3142'),)	{'gene': "B"}
'phiX174'	'CDS'	(('3142', '3312'),)	{'gene': "K"}

$ genbank.py tests/phiX174.gbk -f fasta
>phiX174
gtgtgaggttataacgccgaagcggtaaaaattttaatttttgccgctgagggg
ttgaccaagcgaagcgcggtaggttttctgcttaggagtttaatcatgtttcag

$ genbank.py tests/phiX174.gbk -f fna
>phiX174_CDS_[100..627] [gene="G"]
atgtttcagacttttatttctcgccataattcaaactttttttctgataag
>phiX174_CDS_[636..1622] [gene="H"]
atgtttggtgctattgctggcggtattgcttctgctcttgctggtggcgcc
>phiX174_CDS_[1659..3227]

$ genbank.py tests/phiX174.gbk -f faa
>phiX174_CDS_[100..627] [gene="G"]
MFQTFISRHNSNFFSDKLVLTSVTPASSAPVLQTPKATSSTLYFDSLTVNA
>phiX174_CDS_[636..1622] [gene="H"]
MFGAIAGGIASALAGGAMSKLFGGGQKAASGGIQGDVLATDNNTVGMGDAG
>phiX174_CDS_[1659..3227] [gene="A"]

$ genbank.py tests/phiX174.gbk -f coverage
phiX174	0.965
```



