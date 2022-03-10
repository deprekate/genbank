import io
import gzip
import tempfile

from genbank.locus import Locus

class File(dict):
	def __init__(self, filename=None):
		''' use tempfiles since using next inside a for loop is easier'''
		temp = tempfile.TemporaryFile()
		
		lib = gzip if filename.endswith(".gz") else io
		with lib.open(filename, mode="rb") as fp:
			for line in fp:
				temp.write(line)
				if line.startswith(b'//'):
					temp.seek(0)
					locus = self.parse_locus(temp)
					self[locus.name] = locus
					temp.seek(0)
					temp.truncate()
		temp.close()

	def features(self, include=None, exclude=None):
		for locus in self.values():
			for feature in locus.features(include=include):
				yield feature
	
	def dna(self):
		dna = ""
		for locus in self.values():
			dna += locus.dna
		return dna

	def construct_locus(self):
		'''this method allows for a Locus class to be modified through inheritance in other code '''
		return Locus()

	def parse_locus(self, fp):
		locus = self.construct_locus()
		in_features = False
		current = None

		for line in fp:
			line = line.decode("utf-8")
			if line.startswith('LOCUS'):
				locus.name = line.split()[1]
			elif line.startswith('ORIGIN'):
				in_features = False
				locus.dna = ''
			elif line.startswith('FEATURES'):
				in_features = True
			elif in_features and not line.startswith(" "):
				key,value = line.split()
				#setattr(self, key, value)
				in_features = False
			elif in_features:
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
			elif locus.dna != False:
				locus.dna += line[10:].rstrip().replace(' ','').lower()
		return locus

	def write(self):
		for name in self:
			self[name].write()

