import io
import sys
import gzip
import tempfile
from itertools import chain as chain

from genbank.locus import Locus

class File(dict):
	def __init__(self, filename=None):
		if not hasattr(self, 'locus'):
			self.locus = Locus
		''' use tempfiles since using next inside a for loop is easier'''
		temp = tempfile.TemporaryFile()
		
		lib = gzip if filename.endswith(".gz") else io
		with lib.open(filename, mode="rb") as fp:
			for line in chain(fp, [b'>']):
				# FASTA
				if line.startswith(b'>') and temp.tell() > 1:
					locus = self.parse_locus(temp)
					self[locus.name()] = locus
				temp.write(line)
				# GENBANK
				if line.startswith(b'//'):
					locus = self.parse_locus(temp)
					self[locus.name()] = locus
		temp.close()
	
	def __init_subclass__(cls, locus, **kwargs):
		'''this method allows for a Locus class to be modified through inheritance in other code '''
		super().__init_subclass__(**kwargs)
		cls.locus = locus

	def __iter__(self):
		# hopefully this doesnt break stuff
		return iter(self.values())

	def features(self, include=None, exclude=None):
		for locus in self.values():
			for feature in locus.features(include=include):
				yield feature
	
	def dna(self):
		dna = ""
		for locus in self.values():
			dna += locus.dna
		return dna

	def parse_locus(self, fp):
		locus = self.locus()
		current = None
		field = None
		fasta = False

		fp.seek(0)
		for line in fp:
			line = line.decode("utf-8")
			if line.startswith('\n'):
				pass
			elif not line.startswith(' '):
				field,*value = line.split(maxsplit=1)
				if line.startswith('>'):
					setattr(locus, 'LOCUS', field[1:])
					fasta = True
				elif fasta:
					locus.dna += line.rstrip().replace(' ','').lower()
				elif line.startswith('//'):
					break
				else:
					setattr(locus, field, ''.join(map(str,value)).rstrip())
			elif field == 'FEATURES':
				line = line.rstrip()
				if not line.startswith(' ' * 21):
					while line.endswith(','):
						line += next(fp).decode('utf-8').strip()
					current = locus.read_feature(line)
				else:
					while line.count('"') == 1:
						line += next(fp).decode('utf-8').strip()
					tag,_,value = line[22:].partition('=')
					current.tags[tag] = value #.replace('"', '')
			elif field == 'ORIGIN':
				locus.dna += line.split(maxsplit=1)[1].rstrip().replace(' ','').lower()
			else:
				if not line.startswith('            '):
					field,*value = line.split(maxsplit=1)
					setattr(locus, field, ''.join(map(str,value)).rstrip())
				else:
					value = getattr(locus,field) + line[11:].rstrip()
					setattr(locus,field,value)
		fp.seek(0)
		if fp.writable():
			fp.truncate()
		return locus

	def write(self, outfile=sys.stdout):
		for locus in self:
			locus.write(outfile)

