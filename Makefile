
all:
	python3 -m pip install ../genbank/ --user

clean:
	rm -fr build/
	rm -fr dist/
	rm -fr genbank.egg-info/
	python3 -m pip uninstall -y genbank

upload:
	python3 -m build
	python3 -m twine upload dist/*

build:
	python3 -m build
	python3 -m twine upload dist/*
