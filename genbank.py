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

sys.path.pop(0)
from genbank.file import File

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
	parser.add_argument('-f', '--format', help='Output the features in the specified format', type=str, default='tabular', choices=['tabular','genbank','fasta', 'fna','faa', 'coverage'])
	parser.add_argument('-s', '--slice', help='', type=str, default=None)
	args = parser.parse_args()

	genbank = File(args.infile)
	
	if args.slice:
		left,right = map(int, args.slice.split('..'))
		#print(genbank.current.seq(left, right))
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
			args.outfile.write( getattr(locus, args.format)() )

