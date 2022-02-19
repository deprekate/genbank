from itertools import zip_longest, chain
import textwrap


def rev_comp(dna):
	a = 'acgtrykmbvdh'
	b = 'tgcayrmkvbhd'
	tab = str.maketrans(a,b)
	return dna.translate(tab)[::-1]

def grouper(iterable, n, fillvalue=None):
	args = [iter(iterable)] * n
	return zip_longest(*args, fillvalue=fillvalue)



class Feature():
	def __init__(self, type_, strand, pairs, locus):
		#super().__init__(locus.locus, locus.dna)
		self.type = type_
		self.strand = strand
		# tuplize the pairs
		self.pairs = tuple([tuple(pair) for pair in pairs])
		self.locus = locus
		self.tags = dict()
		self.dna = ''
		self.partial = False

	def frame(self, end):
		if self.type != 'CDS':
			return 0
		elif end == 'right':
			return ((self.right()%3)+1) * self.strand
		elif end == 'left':
			return ((self.left()%3)+1) * self.strand

	def hypothetical(self):
		function = self.tags['product'] if 'product' in self.tags else ''
		if 'hypot'  in function or \
		   'etical' in function or \
		   'unchar' in function or \
		   ('orf' in function and 'orfb' not in function):
			return True
		else:
			return False

	def left(self):
		return int(self.pairs[0][0])
	
	def right(self):
		return int(self.pairs[-1][-1])

	def __str__(self):
		"""Compute the string representation of the feature."""
		return "%s\t%s\t%s\t%s" % (
				repr(self.locus.locus),
				repr(self.type),
				repr(self.pairs),
				repr(self.tags))

	def __repr__(self):
		"""Compute the string representation of the feature."""
		return "%s(%s, %s, %s, %s)" % (
				self.__class__.__name__,
				repr(self.locus),
				repr(self.type),
				repr(self.pairs),
				repr(self.tags))
	def __hash__(self):
		return hash(self.pairs)
	#def __eq__(self, other):
	#	return self.pairs == other.pairs()

	def __lt__(self, other):
		return self.left() < other.left()

	def base_locations(self, full=False):
		if full and self.partial == 'left': 
			for i in range(-((3 - len(self.dna) % 3) % 3), 0, 1):
				yield i+1
		for left,right in self.pairs:
			#left,right = map(int, [ item.replace('<','').replace('>','') for item in self.pair ] )
			for i in range(left,right+1):
				yield i

	def codon_locations(self):
		assert self.type == 'CDS'
		for triplet in grouper(self.base_locations(full=True), 3):
			if triplet[0] >= 1:
				yield triplet

	def codons(self):
		assert self.type == 'CDS'
		for locations in self.codon_locations():
			if self.strand > 0:
				yield ''.join([self.locus.dna[loc-1] if loc else '' for loc in locations])
			else:
				yield rev_comp(''.join([self.locus.dna[loc-1] if loc else '' for loc in locations]))
	
	def translation(self):
		global translate
		aa = []
		codon = ''
		first = 0 if not self.partial else len(self.dna) % 3
		for i in range(first, len(self.dna), 3):
			codon = self.dna[ i : i+3 ]
			if self.strand > 0:
				aa.append(translate.codon(codon))
			else:
				aa.append(translate.codon(rev_comp(codon)))
		if self.strand < 0:
			aa = aa[::-1]
		if aa[-1] in '#*+':
			aa.pop()
		#aa[0] = 'M'
		return "".join(aa)

	def write(self, outfile):
		outfile.write('     ')
		outfile.write( self.type.ljust(16) )
		if not self.strand > 0:
			outfile.write('complement(')
		# the pairs
		if len(self.pairs) > 1:
			outfile.write('join(')
		pairs = []
		for left, right in self.pairs:
			left = max(1,left)
			pair = str(left) + '..' + str(right+2)
			pairs.append(pair)
		outfile.write(','.join(pairs))
		if len(self.pairs) > 1:
			outfile.write(')')
		# the pairs
		if not self.strand > 0:
			outfile.write(')')
		outfile.write('\n')
		for key,value in self.tags.items():
			for line in textwrap.wrap( '/' + str(key) + '=' + str(value) , 58):
				outfile.write('                     ')
				outfile.write(line)
				outfile.write('\n')

