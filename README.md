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
# or add manually by specifing the type,strand,location
locus.add_feature('CDS',+1,[['10','20']])
locus.write()
```
which gives:
```
LOCUS       test                      20 bp
FEATURES             Location/Qualifiers
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
phiX174	1.1168
```
---
You can also *slice* the locus to a specified range and only the nucleotides and
features that occur within the slice are kept. The command to take the first two
hundred bases of the phiX174 genome is shown below. 
```bash
$ genbank.py tests/phiX174.gbk --slice 1..200
LOCUS       phiX174                  200 bp    DNA             PHG
DEFINITION  phiX174
FEATURES             Location/Qualifiers
     rep_origin      13..56
     CDS             100..>200
                     /gene="G"
                     /note="merged"
ORIGIN
        1 gtgtgaggtt ataacgccga agcggtaaaa attttaattt ttgccgctga ggggttgacc
       61 aagcgaagcg cggtaggttt tctgcttagg agtttaatca tgtttcagac ttttatttct
      121 cgccataatt caaacttttt ttctgataag ctggttctca cttctgttac tccagcttct
      181 tcggcacctg ttttacagac
//
```
The *slice* arguement can be paired with all the other output format options:
```bash
$ genbank.py tests/phiX174.gbk --slice 1..200 -f coverage
0.51
```
---
You can also easily edit features by running multiple commands.

Print out the features of a specified **key**:**tag** into a file
```
$ genbank.py tests/phiX174.gbk -k CDS:gene > labels.tsv
```
Change the H of the second gene to something more informative:
(ideally you will have columns from other sources, like excel)
```
perl -pi -e 's/H/Minor spike/' labels.tsv
```

Edit all the features of a specified **key**:**tag** with the 
updated labels:
```
$ genbank.py tests/phiX174.gbk -e CDS:gene < labels.tsv | head
LOCUS       phiX174                 5386 bp    DNA      PHG
FEATURES             Location/Qualifiers
     source          1..5386
     rep_origin      13..56
     CDS             100..627
                     /gene="G"
     CDS             636..1622
                     /gene="Minor spike"
     CDS             1659..3227
```


