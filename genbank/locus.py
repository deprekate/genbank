import re
import sys
import textwrap

from genbank.codons import Last
from genbank.codons import Next
from genbank.codons import Codons
from genbank.feature import Feature
from genbank.translate import Translate

def rev_comp(dna):
	a = 'acgtrykmbvdh'
	b = 'tgcayrmkvbhd'
	tab = str.maketrans(a,b)
	return dna.translate(tab)[::-1]


class Locus(dict):
	def __init__(self, name=None, dna=''):
		if not hasattr(self, 'feature'):
			self.feature = Feature
		self.name = name
		self.dna = dna.lower()
		#self.codons = dict()
		self.translate = Translate()

	def __init_subclass__(cls, feature=Feature, **kwargs):
		'''this method allows for a Feature class to be modified through inheritance in other code '''
		super().__init_subclass__(**kwargs)
		cls.feature = feature

	def fasta(self):
		return ">" + self.name + "\n" + self.seq() + "\n"

	def seq(self, left=1, right=None, strand=1):
		if right is None:
			right = self.length()
		if strand > 0:
			return self.dna[left-1 : right]
		else:
			return rev_comp(self.dna[left-1 : right])

	def length(self):
		return len(self.dna)

	def gc_content(self, seq=None):
		if seq is not None:
			#a = seq.count('a') + seq.count('A')
			c = seq.count('c') + seq.count('C')
			g = seq.count('g') + seq.count('G')
			#t = seq.count('t') + seq.count('T')
			return (c+g) / len(seq) #(a+c+g+t)
		elif not hasattr(self, "gc"):
			#a = self.dna.count('a') + self.dna.count('A')
			c = self.dna.count('c') + self.dna.count('C')
			g = self.dna.count('g') + self.dna.count('G')
			#t = self.dna.count('t') + self.dna.count('T')
			self.gc = (c+g) / len(self.dna) # (a+c+g+t)
		return self.gc

	def pcodon(self, codon):
		codon = codon.lower()
		seq = self.dna + rev_comp(self.dna)
		p = dict()
		p['a'] = seq.count('a') / len(seq)
		p['c'] = seq.count('c') / len(seq)
		p['g'] = seq.count('g') / len(seq)
		p['t'] = seq.count('t') / len(seq)
		return p[codon[0]] * p[codon[1]] * p[codon[2]]

	def rbs(self):
		for feature in self:
			if feature.type == 'CDS':
				if feature.strand > 0:
					start = feature.left()+2
					feature.tags['rbs'] = self.seq(start-30,start)
				else:
					start = feature.right()
					feature.tags['rbs'] = rev_comp(self.seq(start,start+30))
	
	def features(self, include=None, exclude=None):
		for feature in self:
			if not include or feature.type in include:
				yield feature

	def add_feature(self, key, strand, pairs):
		"""Add a feature to the factory."""
		feature = self.feature
		feature = feature(key, strand, pairs, self)
		if feature not in self:
			self[feature] = len(self)
		return feature

	def read_feature(self, line):
		"""Add a feature to the factory."""
		key = line.split()[0]
		#partial  = 'left' if '<' in line else ('right' if '>' in line else False)
		strand = -1 if 'complement' in line else 1
		# this is for weird malformed features
		if ',1)' in line:
			line = line.replace( ",1)" , ",1..1)" )
		#pairs = [pair.split('..') for pair in re.findall(r"<*\d+\.\.>*\d+", line)]
		#pairs = [map(int, pair.split('..')) for pair in re.findall(r"<?\d+\.{0,2}>?\d+", line.replace('<','').replace('>','') )]
		pairs = [ pair.split('..') for pair in re.findall(r"<?\d+\.{0,2}>?\d*", line) ]
		# tuplize the pairs
		pairs = tuple([tuple(pair) for pair in pairs])
		feature = self.add_feature(key, strand, pairs)
		return feature

	def gene_coverage(self):
		''' This calculates the protein coding gene coverage, which should be around 1 '''
		cbases = tbases = 0	
		for locus in self.values():
			dna = [False] * len(locus.dna)
			seen = dict()
			for feature in locus.features(include=['CDS']):
				for i in feature.codon_locations():
					dna[i-1] = True
			cbases += sum(dna)
			tbases += len(dna)
		return 3 * cbases / tbases

	def write(self, outfile=sys.stdout):
		outfile.write('LOCUS	   ')
		outfile.write(self.name)
		outfile.write(str(len(self.dna)).rjust(10))
		outfile.write(' bp	DNA			 UNK')
		outfile.write('\n')
		outfile.write('DEFINITION  ' + self.name + '\n')
		outfile.write('FEATURES			 Location/Qualifiers\n')
		#outfile.write('	 source		  1..')
		#outfile.write(str(len(self.dna)))
		#outfile.write('\n')
		for feature in self:
			feature.write(outfile)
		outfile.write('ORIGIN')
		i = 0
		dna = textwrap.wrap(self.dna, 10)
		for block in dna:
			if(i%60 == 0):
				outfile.write('\n')
				outfile.write(str(i+1).rjust(9))
				outfile.write(' ')
				outfile.write(block.lower())
			else:
				outfile.write(' ')
				outfile.write(block.lower())
			i += 10
		outfile.write('\n')
		outfile.write('//')
		outfile.write('\n')

	def last(self, n, strand, codons):
		if isinstance(codons, str):
			codons = [codons]
		if strand < +1:
			codons = list(map(rev_comp, codons))
		for i in range(n, 0, -3):
			if self.dna[i:i+3] in codons:
				return i
		return None

	def next(self, n, strand, codons):
		if isinstance(codons, str):
			codons = [codons]
		if strand < +1:
			codons = list(map(rev_comp, codons))
		for i in range(n, self.length(), 3):
			#print(list(codons), self.dna[i:i+3])
			if self.dna[i:i+3] in codons:
				return i
		return None

	def nearest(self, n, strand, codons):
		_last = self.last(n,strand,codons)
		_next = self.next(n,strand,codons)
		if n - _last < _next - n:
			return _last
		else:
			return _next

	def distance(self, n, strand, codons):
		nearest = self.nearest(n, strand, codons)
		return n - nearest if nearest < n else nearest - n

	def codon_rarity(self, codon):
		if not hasattr(self, 'rarity'):
			self.rarity = {a+b+c : 1 for a in 'acgt' for b in 'acgt' for c in 'acgt'}
			for feature in self:
				if feature.type == 'CDS':
					for _codon in feature.codons():
						if _codon in self.rarity:
							self.rarity[_codon] += 1
		total = sum(self.rarity.values())
		self.rarity = {codon:self.rarity[codon]/total for codon in self.rarity}
		return round(self.rarity[codon], 5)

		



		

