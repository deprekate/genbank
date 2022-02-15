
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

from genbank.locus import Locus


def nint(x):
    return int(x.replace('<','').replace('>',''))


class Read(Locus):
	''' inherit a locus class to fill '''
	def __init__(self, fp):
		self.locus = None
		self.dna = False
		in_features = False
		current = None

		for line in fp:
			line = line.decode("utf-8")
			if line.startswith('LOCUS'):
				self.locus = line.split()[1]
			elif line.startswith('ORIGIN'):
				in_features = False
				self.dna = ''
			elif line.startswith('FEATURES'):
				in_features = True
			elif in_features:
				line = line.rstrip()
				if not line.startswith(' ' * 21):
					while line.endswith(','):
						line += next(fp).decode('utf-8').strip()
					current = self.read_feature(line)
				else:
					while line.count('"') == 1:
						line += next(fp).decode('utf-8').strip()
					tag,_,value = line[22:].partition('=')
					current.tags[tag] = value.replace('"', '')
			elif self.dna != False:
				self.dna += line[10:].rstrip().replace(' ','').lower()

	def features(self, include=None, exclude=None):
		for  feature in self:
			if not include or feature.type in include:
				yield feature

	def read_feature(self, line):
		"""Add a feature to the factory."""
		key = line.split()[0]
		partial  = 'left' if '<' in line else ('right' if '>' in line else False)
		strand = -1 if 'complement' in line else 1
		pairs = [pair.split('..') for pair in re.findall(r"<*\d+\.\.>*\d+", line)]
		# this is for weird malformed features
		if ',1)' in line:
			pairs.append(['1','1'])
		# tuplize the pairs
		pairs = tuple([tuple(pair) for pair in pairs])
		feature = self.add_feature(key, strand, pairs)
		feature.partial = partial
		return feature
		














