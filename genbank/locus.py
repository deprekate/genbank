import sys

from genbank.feature import Feature

def rev_comp(dna):
	a = 'acgtrykmbvdh'
	b = 'tgcayrmkvbhd'
	tab = str.maketrans(a,b)
	return dna.translate(tab)[::-1]


class Locus(dict):
	def __init__(self, locus, dna):
		self.locus = locus
		self.dna = dna.lower()
		self.locations = cd.Locations(self.dna)
		
		seq = self.dna + rev_comp(self.dna)
		length = len(seq)
		self.p = dict()
		self.p['a'] = seq.count('a') / length
		self.p['c'] = seq.count('c') / length
		self.p['g'] = seq.count('g') / length
		self.p['t'] = seq.count('t') / length

		self.wobble = {'gcc':'cgg', 'gct':'cgg', 'gca':'cgt', 'gcg':'cgt',
					   'aga':'tct', 'agg':'tct', 'cga':'gct', 'cgg':'gct', 'cgt':'gcg', 'cgc':'gcg',
					   'gac':'ctg', 'gat':'ctg', 
					   'aac':'ttg', 'aat':'ttg',
					   'tgc':'acg', 'tgt':'acg',
					   'gaa':'ctt', 'gag':'ctt',
					   'caa':'gtt', 'cag':'gtt',
					   'gga':'cct', 'ggg':'cct', 'ggc':'ccg', 'ggt':'ccg',
					   'cac':'gtg', 'cat':'gtg',
					   'ata':'tat', 'atc':'tag', 'att': 'tag',
					   'tta':'aat', 'ttg':'aat', 'cta':'gta', 'ctg':'gta', 'ctt':'gtg', 'ctc':'gtg',
					   'aaa':'ttt', 'aag':'ttt',
					   'atg':'tac',
					   'ttc':'aag', 'ttt':'aag',
					   'cca':'ggt', 'ccg':'ggt', 'cct':'ggg', 'ccc':'ggg',
					   'agc':'tcg', 'agt':'tcg', 'tca':'atg', 'tcg':'atg', 'tcc':'agg', 'tct':'agg',
					   'aca':'tgt', 'acg':'tgt', 'acc':'tgg', 'act':'tgg',
					   'tgg':'acc',
					   'tac':'atg', 'tat':'atg',
					   'gta':'cat', 'gtg':'cat', 'gtc':'cag', 'gtt':'cag',
					   'tag':'atc', 'taa':'att', 'tga':'act'
					   }

	def seq(self, left, right):
		return self.dna[left-1 : right]

	def length(self):
		return len(self.dna)

	def pcodon(self, codon):
		codon = codon.lower()
		return self.p[codon[0]] * self.p[codon[1]] * self.p[codon[2]]

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
		feature = Feature(key, strand, pairs, self)
		if feature not in self:
			self[feature] = True
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

	def write(self, outfile):
		outfile.write('LOCUS       ')
		outfile.write(self.locus)
		outfile.write(str(len(self.dna)).rjust(10))
		outfile.write(' bp    DNA             UNK')
		outfile.write('\n')
		outfile.write('DEFINITION  ' + self.locus + '\n')
		outfile.write('FEATURES             Location/Qualifiers\n')
		outfile.write('     source          1..')
		outfile.write(str(len(self.dna)))
		outfile.write('\n')

		for feature in sorted(self):
			feature.write(outfile)
		outfile.write('//')
		outfile.write('\n')


			
