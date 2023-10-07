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
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=argparse.RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file in genbank format')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-f', '--format', help='Output the features in the specified format', type=str, default='genbank', choices=File.formats)
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
		perfect = partial = total = fp = 0
		compare = File(args.compare)
		partial = dict()
		perfect = dict()
		mistake = dict()
		for locus,other in zip(genbank,compare):
			pairs = dict()
			for feature in locus.features(include='CDS'):
				if feature.strand > 0:
					pairs[feature.pairs[-1][-1]] = feature.pairs[ 0][ 0]
				else:
					pairs[feature.pairs[ 0][ 0]] = feature.pairs[-1][-1]
			total += len(pairs)
			for feature in other.features(include='CDS'):
				start,stop = feature.pairs[ 0][ 0],feature.pairs[-1][-1]
				if feature.strand < 0:
					stop,start = start,stop
				if stop in pairs:
					partial[stop] = True
					if start == pairs[stop]:
						perfect[stop] = True
				else:
					mistake[stop] = True
					#print(feature)
		partial = len(partial)
		perfect = len(perfect)
		mistake = len(mistake)
		args.outfile.write(f"{partial}\t{partial/total}\t{perfect}\t{perfect/total}\t{total}\t{mistake}\n")
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
				elif args.add == 'genemarks':
					if line.startswith(' ') and 'Gene' not in line and '#' not in line:
						n,strand,left,right,*_ = line.split()
						locus.add_feature('CDS',strand,[[left,right]],{'note':['genemarkS']})
				elif args.add == 'genemark':
					if line.startswith(' ') and (' direct ' in line or ' complement ' in line):
						left,right,strand,*_ = line.split()
						strand = '+' if ('direct' in line or '+' in line) else '-'
						locus.add_feature('CDS',strand,[[left,right]],{'note':['genemark']})
					elif line.startswith('List of Regions of interest'):
						break
				elif args.add == 'glimmer':
					if not line.startswith('>'):
						_,left,right,(strand,*_),*_ = line.split()
						if strand == '-':
							left, right = right, left
						if int(left) > int(right):
							if strand == '-':
								pairs = [['1', right], [left, str(locus.length())]]
							else:
								pairs = [[left, str(locus.length())], ['1', right]]
						else:
							pairs = [[left,right]]
						locus.add_feature('CDS',strand,pairs,{'note':['glimmer3']})
				elif args.add == 'gff':
					if not line.startswith('#') and len(line) > 2:
						try:
							_,_,key,left,right,_,strand,_,tags,*_ = line.rstrip('\n').split('\t')
						except:
							print('Error:')
							print(line)
							exit()
						key = key[:14]
						tags = {key: ['"%s"' % '='.join(val)] for tag in tags.split(';') for key,*val in [tag.split('=')]}
						locus.add_feature(key,strand,[[left,right]],tags)
				else:
					if not line.startswith('#') and len(line) > 2:
						mapping = dict()
						for key,val in zip(args.add.split(','),line.split()):
							mapping[key] = val
						left,right,strand = [mapping[key] for key in ['left','right','strand']]
						locus.add_feature('CDS',strand,[[left,right]],{})

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
			left -= 1
		elif ':' in args.slice:
			left,right = args.slice.split(':')
			if '+' in right and '-' in right:
				left = right = eval(left + right)
			elif '+' in right:
				right = eval(left + right)
			elif '-' in right:
				left,right = eval(left + right) , left
			left,right = map(int, [left,right])
		elif '-' in args.slice:
			left,right = map(int, args.slice.split('-'))
			right += 1
		else:
			raise Exception("re-circularization not implemented yet")
			left = int(args.slice)
			right = left+1
		for locus in genbank:
			locus = locus.slice(left,right)
	if args.key:
		key,qualifier = args.key.replace('/',':').split(':')
		for feature in genbank.features(include=key):
			args.outfile.print('\t'.join(feature.tags.get(qualifier,'')))
			args.outfile.print("\n")
	elif args.format == 'genbank':
		if args.revcomp:
			raise Exception("not implemented yet")
		genbank.write(args)	
	elif args.format in ['fna','faa','gff', 'gff3','tabular']:
		for locus in genbank:
			locus.write(args)
	elif args.format in ['fasta']:
		# ONCE THE REVCOMP IS IMPLEMENTED FOR GENBANK MOST OF THESE
		# FORMATS WILL CONDENSE INTO A SINGLE BLOCK
		for locus in genbank:
			if args.revcomp:
				locus.dna = locus.seq(strand=-1)
			locus.write(args)
	elif args.format == 'bases':
		strand = -1 if args.revcomp else +1
		for locus in genbank:
			args.outfile.print(locus.seq(strand=strand))
			args.outfile.print('\n')
	elif args.format == 'rarity':
		rarity = dict()
		for locus in genbank:
			for codon,freq in sorted(locus.codon_rarity().items(), key=lambda item: item[1]):
				args.outfile.print(codon)
				args.outfile.print('\t')
				args.outfile.print(round(freq,5))
				args.outfile.print('\n')
	elif args.format == 'coverage':
		cbases = tbases = 0
		for locus in genbank:
			c,t = locus.gene_coverage()
			cbases += c
			tbases += t
		#args.outfile.print( name )
		#args.outfile.print( '\t' )
		args.outfile.print( cbases / tbases )
		args.outfile.print( '\n' )
	elif args.format in ['gc','gcfp']:
		for locus in genbank:
			args.outfile.print(locus.name())
			args.outfile.print('\t')
			if args.format == 'gc':
				args.outfile.print(locus.gc_content())
			else:
				args.outfile.print(locus.gc_fp())
			args.outfile.print('\n')
	elif args.format == 'taxonomy':
		for locus in genbank:
			if 'SOURCE' in locus.groups:
				s = locus.groups['SOURCE'][0].replace('\n','\t').replace(' '*12,'').replace(';\t','; ')
				args.outfile.print(s)
				args.outfile.print('\n')
	elif args.format in ['split']:
		folder = args.outfile.name if args.outfile.name != '<stdout>' else ''
		for locus in genbank:
			with open(os.path.join(folder,name + '.fna'), 'w') as f:
				f.write('>')
				f.write(locus.name())
				f.write('\n')
				f.write(locus.seq())
				f.write('\n')
	elif args.format == 'testcode':
		for locus in genbank:
			args.outfile.print(locus.name())
			args.outfile.print('\t')
			args.outfile.print(locus.testcode())
			args.outfile.print('\n')




