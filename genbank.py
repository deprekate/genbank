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
import fileinput

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

def _print(self, item):
    if isinstance(item, str):
        self.write(item)
    else:
        self.write(str(item))

if __name__ == "__main__":
	choices = 	['tabular','genbank','fasta', 'fna','faa', 'coverage','rarity','bases','gc','gcfp', 'taxonomy','part', 'gff', 'gff3', 'testcode']
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file in genbank format')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-f', '--format', help='Output the features in the specified format', type=str, default='genbank', choices=choices)
	parser.add_argument('-s', '--slice', help='This slices the infile at the specified coordinates. \nThe range can be in one of three different formats:\n    -s 0-99      (zero based string indexing)\n    -s 1..100    (one based GenBank indexing)\n    -s 50:+10    (an index and size of slice)', type=str, default=None)
	parser.add_argument('-g', '--get', action="store_true")
	parser.add_argument('-r', '--revcomp', action="store_true")
	parser.add_argument('-a', '--add', help='This adds features the shell input via < features.txt', type=str, default=None)
	parser.add_argument('-e', '--edit', help='This edits the given feature key with the value from the shell input via < new_keys.txt', type=str, default=None)
	parser.add_argument('-k', '--key', help='Print the given keys [and qualifiers]', type=str, default=None)
	parser.add_argument('-c', '--compare', help='Compares the CDS of two genbank files', type=str, default=None)
	args = parser.parse_args()
	args.outfile.print = _print.__get__(args.outfile)

	if not args.get:
		genbank = File(args.infile)
	else:
		#raise Exception("not implemented yet")
		# not ready yet
		accession,rettype = args.infile.split('.')
		url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + accession + '&rettype=' + rettype + '&retmode=text'
		with urllib.request.urlopen(url) as response:
			with tempfile.NamedTemporaryFile() as tmp:
				shutil.copyfileobj(response, tmp)
				genbank = File(tmp.name)
		
	if args.compare:
		perfect = partial = total = 0
		compare = File(args.compare)
		for (name,locus),( _ ,other) in zip(genbank.items(),compare.items()):
			pairs = dict()
			for feature in locus.features(include='CDS'):
				if feature.strand > 0:
					pairs[feature.pairs[-1][-1]] = feature.pairs[ 0][ 0]
				else:
					pairs[feature.pairs[ 0][ 0]] = feature.pairs[-1][-1]
			total += len(pairs)
			for feature in other.features(include='CDS'):
				if feature.strand > 0:
					if feature.pairs[-1][-1] in pairs:
						partial += 1
						if feature.pairs[ 0][ 0] == pairs[feature.pairs[-1][-1]]:
							perfect += 1
				else:
					if feature.pairs[ 0][ 0] in pairs:
						partial += 1
						if feature.pairs[-1][-1] == pairs[feature.pairs[ 0][ 0]]:
							perfect += 1
		args.outfile.print(partial)
		args.outfile.print('\t')
		args.outfile.print('(')
		args.outfile.print(partial/total)
		args.outfile.print(')')
		args.outfile.print('\t')
		args.outfile.print(perfect)
		args.outfile.print('\t')
		args.outfile.print('(')
		args.outfile.print(perfect/total)
		args.outfile.print(')')
		args.outfile.print('\t')
		args.outfile.print(total)
		args.outfile.print('\n')
		exit()
	if args.add:
		# this only works for single sequence files
		stdin = []
		if not sys.stdin.isatty():
			stdin = sys.stdin.readlines()
		for locus in genbank:
			for line in stdin:
				if args.add == 'genbank':
					pass
				elif args.add == 'genemark':
					if line.startswith(' ') and 'Gene' not in line and '#' not in line:
						key = 'CDS'
						n,strand,left,right,*_ = line.split()
						locus.add_feature(key,strand,[[left,right]],{'note':['genemarkS']})
				elif args.add == 'glimmer':
					if not line.startswith('>'):
						key = 'CDS'
						n,left,right,(strand,*_),*_ = line.split()
						if strand == '-':
							left, right = right, left
							pairs = [[left, right]] if int(left) > int(right) else [[left, str(locus.length())], ['1', right]]
						else:
							pairs = [[left, right]] if int(left) > int(right) else [[left, right]]
						locus.add_feature(key,strand,pairs,{'note':['glimmer3']})
				elif args.add == 'gff':
					if not line.startswith('#'):
						name,other,key,left,right,_,strand,_,tags = line.rstrip('\n').split('\t')
						tags = dict((key,value) for tag in tags.split(';') for key,*value in [tag.split('=')])
						locus.add_feature(key,strand,[[left,right]],tags)
	elif args.edit:
		if not sys.stdin.isatty():
			stdin = sys.stdin.readlines()
			#sys.stdin = open('/dev/tty')
		key,qualifier = args.edit.replace('/',':').split(':')
		for feature,values in zip(genbank.features(include=[key]), stdin):
			feature.tags[qualifier] = list()
			for value in values.rstrip().split('\t'):
				feature.tags[qualifier].append(value)
	if args.slice:
		if '..' in args.slice:
			left,right = map(int, args.slice.split('..'))
			left = left - 1
		elif ':' in args.slice:
			left,right = args.slice.split(':')
			if '+' in right and '-' in right:
				left = eval(left + right)
				right = eval(left + right)
			elif '+' in right:
				right = eval(left + right)
			elif '-' in right:
				left,right = eval(left + right) , left
			left,right = map(int, [left,right])
		elif '-' in args.slice:
			left,right = map(int, args.slice.split('-'))
			right = right + 1
		else:
			raise Exception("re-circularization not implemented yet")
			left = int(args.slice)
			right = left+1
		for name,locus in genbank.items():
			locus = locus.slice(left,right)
	if args.key:
		key,qualifier = args.key.replace('/',':').split(':')
		for feature in genbank.features(include=key):
			args.outfile.print('\t'.join(feature.tags[qualifier]))
			args.outfile.print("\n")
	elif args.format == 'genbank':
		if args.revcomp:
			raise Exception("not implemented yet")
		genbank.write(args.outfile)	
	elif args.format == 'tabular':
		for feature in genbank.features(include=['CDS']):
			args.outfile.print(feature)
			args.outfile.print("\t")
			args.outfile.print(feature.seq())
			args.outfile.print("\n")
	elif args.format in ['gff', 'gff3']:
		for locus in genbank:
			locus.write(args.outfile, args)
	elif args.format in ['fna','faa']:
		for name,locus in genbank.items():
			for feature in locus.features(include=['CDS']):
				args.outfile.print( getattr(feature, args.format)() )
	elif args.format in ['fasta']:
		for name,locus in genbank.items():
			if args.revcomp:
				locus.dna = locus.seq(strand=-1)
			args.outfile.print( getattr(locus, args.format)() )
	elif args.format == 'coverage':
		cbases = tbases = 0
		for name,locus in genbank.items():
			c,t = locus.gene_coverage()
			cbases += c
			tbases += t
		#args.outfile.print( name )
		#args.outfile.print( '\t' )
		args.outfile.print( cbases / tbases )
		args.outfile.print( '\n' )
	elif args.format == 'rarity':
		rarity = dict()
		for name,locus in genbank.items():
			for codon,freq in sorted(locus.codon_rarity().items(), key=lambda item: item[1]):
				args.outfile.print(codon)
				args.outfile.print('\t')
				args.outfile.print(round(freq,5))
				args.outfile.print('\n')
	elif args.format == 'bases':
		strand = -1 if args.revcomp else +1
		for name,locus in genbank.items():
			args.outfile.print(locus.seq(strand=strand))
			args.outfile.print('\n')
	elif args.format in ['gc','gcfp']:
		for name,locus in genbank.items():
			args.outfile.print(locus.name())
			args.outfile.print('\t')
			if args.format == 'gc':
				args.outfile.print(locus.gc_content())
			else:
				args.outfile.print(locus.gc_fp())
			args.outfile.print('\n')
	elif args.format == 'taxonomy':
		for name,locus in genbank.items():
			args.outfile.print(locus.groups['SOURCE'][0].replace('\n','\t').replace('            ','').replace(';\t','; ') )
			args.outfile.print('\n')
	elif args.format in ['part']:
		folder = args.outfile.name if args.outfile.name != '<stdout>' else ''
		for name,locus in genbank.items():
			with open(os.path.join(folder,name + '.fna'), 'w') as f:
				f.write('>')
				f.write(name)
				f.write('\n')
				f.write(locus.seq())
				f.write('\n')
	elif args.format == 'testcode':
		for name,locus in genbank.items():
			args.outfile.print(locus.name())
			args.outfile.print('\t')
			args.outfile.print(locus.testcode())
			args.outfile.print('\n')




