from itertools import zip_longest, chain
import textwrap
import copy

from genbank.sequence import Seq

def rev_comp(dna):
	a = 'acgtrykmbvdh'
	b = 'tgcayrmkvbhd'
	tab = str.maketrans(a,b)
	return dna.translate(tab)[::-1]

def grouper(iterable, n, fillvalue=None):
	args = [iter(iterable)] * n
	return zip_longest(*args, fillvalue=fillvalue)

def nint(s):
	return int(s.replace('<','').replace('>',''))

class Feature():
	def __init__(self, type_, strand, pairs, locus, tags=None):
		#super().__init__(locus.locus, locus.dna)
		self.type = type_
		self.strand = strand
		# tuplize the pairs
		self.pairs = tuple([tuple(pair) for pair in pairs])
		self.locus = locus
		self.tags = tags if tags else dict()
		#self.dna = ''
		#self.partial = False

	def frame(self, end='left'):
		if self.type != 'CDS':
			return 0
		elif end == 'right' or self.partial() == 'left':
			return (self.right()%3+1) * self.strand
		elif end == 'left' or self.partial() == 'right':
			return (self.left()%3+1) * self.strand

	def left(self):
		# convert genbank 1-based indexing to standard 0-based
		# this is probably not right for <1 that should -2
		return nint(self.pairs[0][0]) - 1
	
	def right(self):
		# convert genbank 1-based indexing to standard 0-based
		return nint(self.pairs[-1][-1]) - 3
	
	def length(self):
		return len(self.seq())

	'''
	def seq(self):
		seq = list()
		for pair in self.pairs:
			left,right = map(nint, pair)
			seq.append(self.locus.seq(left-1,right, self.strand))
		if self.strand > 0:
			return ''.join(seq)
		else:
			return ''.join(seq[::-1])
	'''

	def seq(self):
		# lol this is chatgpts more optimized version
		seq = [self.locus.seq(nint(pair[0]) - 1, nint(pair[-1]), self.strand) for pair in self.pairs]
		return ''.join(seq) if self.strand > 0 else ''.join(seq[::-1])

	'''
	def loc(self):
		loc = list()
		for pair in self.pairs:
			left,right = map(nint, pair)
			loc.extend(list(range(left-1,right)))
		return loc
	'''
	def loc(self):
		# lol this is chatgpts more optimized version
		return [val for pair in self.pairs for val in range(nint(pair[0])-1, nint(pair[-1]))]
	
	def codons(self, loc=False):
		assert self.type == 'CDS'
		dna = self.seq()
		# should I return partial codons?
		partial_type = self.partial()
		remainder = len(dna) % 3
		if partial_type == 'left':
			dna = dna[remainder:] if self.strand > 0 else dna[:-remainder]
		elif partial_type == 'right':
			dna = dna[:-remainder] if self.strand > 0 else dna[remainder:]
		if len(dna) % 3:
			# not sure how to handle bad features
			dna = ''
		for triplet in grouper(dna, 3):
			yield ''.join(triplet)
		return

	def codon_locations(self):
		assert self.type == 'CDS'
		loc = self.loc()
		if self.strand < 0:
			loc.reverse()
		# should I return partial codons?
		partial_type = self.partial()
		remainder = len(loc) % 3
		if partial_type == 'left':
			loc = loc[remainder:] if self.strand > 0 else loc[:-remainder]
		elif partial_type == 'right':
			loc = loc[:-remainder] if self.strand > 0 else loc[remainder:]
		if len(loc) % 3:
			# not sure how to handle bad features
			loc = []
		for locs in grouper(loc, 3):
			yield locs
		return

	def fna(self):
		return self.header() + self.seq() + "\n"

	def faa(self):
		return self.header() + self.translation() + "\n"
	
	def header(self):
		header = ">" + self.locus.name() + "_CDS_[" + self.locations() + "]"
		for tag,values in self.tags.items():
			if tag != 'translation':
				for value in values:
					if value:
						header += " [" + tag + "=" + str(value) +"]"
					else:
						header += " [" + tag +"]"
		return header + "\n"

	def hypothetical(self):
		function = ''.join(self.tags['product']).lower() if 'product' in self.tags else ''
		if 'hypot'  in function or \
		   'etical' in function or \
		   'unchar' in function or \
		   ('orf' in function and 'orfb' not in function):
			return True
		else:
			return False

	def partial(self):
		partial_left  = any(['<' in item for pair in self.pairs for item in pair])
		partial_right = any(['>' in item for pair in self.pairs for item in pair])
		if partial_left and partial_right:
			# this really shouldnt even happen, maybe raise an error?
			return 'both'
		elif partial_left:
			return 'left'
		elif partial_right:
			return 'right'
		else:
			return None

	def is_type(self, _type):
		if self.type == _type:
			return True
		else:
			return False

	def is_joined(self):
		if len(self.pairs) > 1:
			return True
		return False

	'''
	def __iter__(self):
		for left,*right in self.pairs:
			if right:
				right = right[0]
			else:
				right = left
			yield nint(left)-1 , nint(right)-1
	'''
	def __str__(self):
		"""Compute the string representation of the feature."""
		return "%s\t%s\t%s\t%s" % (
				repr(self.locus.name()),
				repr(self.type),
				repr(self.pairs),
				repr(self.tags))

	def __repr__(self):
		"""Compute the string representation of the feature."""
		return "%s(%s, %s, %s, %s)" % (
				self.__class__.__name__,
				repr(self.locus.name()),
				repr(self.type),
				repr(self.pairs),
				repr(self.tags))
	def __hash__(self):
		return hash(self.pairs)
	#def __eq__(self, other):
	#	return self.pairs == other.pairs()

	def __lt__(self, other):
		if self.left() == other.left():
			return self.right() < other.right()
		else:
			return self.left() < other.left()

	def locations(self):
		pairs = []
		#for left, *right in self.pairs:
		for pair in self.pairs:
			pairs.append("..".join(pair))
		location = ','.join(pairs)
		if len(pairs) > 1:
			location = 'join(' + location + ')'
		if self.strand < 0:
			location = 'complement(' + location + ')'
		return location


	def split(self):
		a = copy.copy(self)
		b = copy.copy(self)
		return a,b

	def translation(self):
		aa = []
		codon = ''
		first = self.length() % 3 if (self.partial() == 'left' and self.strand > 0) or (self.partial() == 'right' and self.strand < 0) else 0
		dna = self.seq()
		for i in range(first, self.length(), 3):
			codon = dna[ i : i+3 ]
			aa.append(self.locus.translate.codon(codon))
		#if self.strand < 0:
		#	aa = aa[::-1]
		# keeping the stop codon character adds 'information' as does which of
		# the stop codons it is. It is the better way to write the fasta
		#if aa[-1] in '#*+':
		#	aa.pop()
		# keeping the first amino acid also adds 'information' to the fasta
		#aa[0] = 'M'
		return "".join(aa)

	def write(self, outfile):
		outfile.write('     ')
		outfile.write( self.type.ljust(16) )
		if self.strand < 0:
			outfile.write('complement(')
		# the pairs
		if len(self.pairs) > 1:
			outfile.write('join(')
		pairs = []
		#for left, *right in self.pairs:
		for pair in self.pairs:
			#pair = left + '..' + str(right[0]) if right else str(left)
			pairs.append("..".join(pair))
		outfile.write(','.join(pairs))
		if len(self.pairs) > 1:
			outfile.write(')')
		# the pairs
		if self.strand < 0:
			outfile.write(')')
		outfile.write('\n')
		for tag,values in self.tags.items():
			for value in values:
				if value is not None:
					for line in textwrap.wrap( '/' + str(tag) + '=' + str(value) , 58):
						outfile.write('                     ')
						outfile.write(line)
						outfile.write('\n')
				else:
					outfile.write('                     ')
					outfile.write('/' + str(tag))
					outfile.write('\n')

