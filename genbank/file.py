import io
import gzip
import tempfile

from genbank.read import Read

class File(dict):
	def __init__(self, filename, current=None):
		''' use tempfiles since using next inside a for loop is easier'''
		temp = tempfile.TemporaryFile()
		
		lib = gzip if filename.endswith(".gz") else io
		with lib.open(filename, mode="rb") as fp:
			for line in fp:
				temp.write(line)
				if line.startswith(b'//'):
					temp.seek(0)
					locus = Read(temp)
					self[locus.locus] = locus
					temp.seek(0)
					temp.truncate()
		temp.close()

	def features(self, include=None, exclude=None):
		for locus in self.values():
			for feature in locus.features(include=include):
				yield feature

