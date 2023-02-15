class Seq(str):
	# this is just to capture negative string indices as zero
	def __getitem__(self, key):
		if isinstance(key, slice):
			if key.stop < 0:
				key = slice(0, 0, key.step)
			elif key.start < 0:
				key = slice(0, key.stop, key.step)
		elif isinstance(key, int) and key >= len(self):
			return ''
		return super().__getitem__(key)
