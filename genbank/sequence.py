def rev_comp(dna):
	a = 'acgtrykmbvdh'
	b = 'tgcayrmkvbhd'
	tab = str.maketrans(a,b)
	return dna.translate(tab)[::-1]


class Seq(str):
	# this is just to capture negative string indices as zero
	def __getitem__(self, key):
		if isinstance(key, slice):
			if key.stop and key.stop < 0:
				key = slice(0, 0, key.step)
			elif key.start and key.start < 0:
				key = slice(0, key.stop, key.step)
		elif isinstance(key, int) and key >= len(self):
			return ''
		return super().__getitem__(key)
