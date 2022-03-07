
all:
	python3 setup.py install --user

clean:
	rm -fr build/
	rm -fr dist/
	rm -fr genbank.egg-info/
	pip uninstall -y genbank
