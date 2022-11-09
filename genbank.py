#!/usr/bin/env python3
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 
import os
import io
import sys
import gzip
import re
import argparse
import tempfile
from collections import Counter
from argparse import RawTextHelpFormatter
from itertools import zip_longest, chain
import shutil
import tempfile
import urllib.request

sys.path.pop(0)
from genbank.file import File

def get(x):
	return True

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def nint(x):
    return int(x.replace('<','').replace('>',''))


if __name__ == "__main__":
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file in genbank format')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-f', '--format', help='Output the features in the specified format', type=str, default='tabular', choices=['tabular','genbank','fasta', 'fna','faa', 'coverage','rarity','bases','gc','taxonomy','revcomp','part'])
	parser.add_argument('-s', '--slice', help='', type=str, default=None)
	parser.add_argument('-g', '--get', action="store_true")
	parser.add_argument('-r', '--revcomp', action="store_true")
	args = parser.parse_args()

	if not args.get:
		genbank = File(args.infile)
	else:
		raise Exception("not implemented yet")
		# not ready yet
		accession,rettype = args.infile.split('.')
		with urllib.request.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + accession + '&rettype=' + rettype + '&retmode=text') as response:
			with tempfile.NamedTemporaryFile() as tmp:
				shutil.copyfileobj(response, tmp)
				genbank = File(tmp.name)

	if args.slice:
		if '..' in args.slice:
			left,right = map(int, args.slice.split('..'))
			left = left-1
		elif '-' in args.slice:
			left,right = map(int, args.slice.split('-'))
			right = right+1
		elif ':' in args.slice:
			left,right = args.slice.split(':')
			if '+' in right:
				right = eval(left + right)
			left,right = map(int, [left,right])
		else:
			left = int(args.slice)
			right = left+1
		for name,locus in genbank.items():
			locus = locus.slice(left,right)
	if args.format == 'genbank':
		genbank.write(args.outfile)	
	elif args.format == 'tabular':
		for feature in genbank.features(include=['CDS']):
			args.outfile.write(str(feature))
			args.outfile.write("\t")
			args.outfile.write(feature.seq())
			args.outfile.write("\n")
	elif args.format in ['fna','faa']:
		for name,locus in genbank.items():
			for feature in locus.features(include=['CDS']):
				args.outfile.write( getattr(feature, args.format)() )
	elif args.format in ['fasta']:
		for name,locus in genbank.items():
			if args.revcomp:
				locus.dna = locus.seq(strand=-1)
			args.outfile.write( getattr(locus, args.format)() )
	elif args.format == 'coverage':
		for name,locus in genbank.items():
			args.outfile.write( name )
			args.outfile.write( '\t' )
			args.outfile.write( str(locus.gene_coverage()) )
			args.outfile.write( '\n' )
	elif args.format == 'rarity':
		rarity = dict()
		for name,locus in genbank.items():
			for codon,freq in sorted(locus.codon_rarity().items(), key=lambda item: item[1]):
				args.outfile.write(codon)
				args.outfile.write('\t')
				args.outfile.write(str(round(freq,5)))
				args.outfile.write('\n')
	elif args.format == 'bases':
		for name,locus in genbank.items():
			args.outfile.write(locus.dna)
			args.outfile.write('\n')
	elif args.format == 'gc':
		for name,locus in genbank.items():
			args.outfile.write(locus.name())
			args.outfile.write('\t')
			args.outfile.write(str(locus.gc_content()))
			args.outfile.write('\n')
	elif args.format == 'taxonomy':
		for name,locus in genbank.items():
			args.outfile.write(locus.name())
			args.outfile.write('\t')
			args.outfile.write(locus.ORGANISM)
			args.outfile.write('\n')
	elif args.format in ['part']:
		folder = args.outfile.name if args.outfile.name != '<stdout>' else ''
		for name,locus in genbank.items():
			with open(os.path.join(folder,name + '.fna'), 'w') as f:
				f.write('>')
				f.write(name)
				f.write('\n')
				f.write(locus.seq())
				f.write('\n')




