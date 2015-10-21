VMD = vmd

test:
	${VMD} -python -dispdev none -e pyvmd/tests/__init__.py -args discover

coverage:
	python-coverage erase
	-rm -r htmlcov
	-COVERAGE=1 ${VMD} -python -dispdev none -e pyvmd/tests/__init__.py -args discover
	python-coverage html -d htmlcov

pylint:
	-PYTHONPATH="/usr/lib/vmd/scripts/python" pylint pyvmd

pep8:
	-pep8 pyvmd --max-line-length=119

pepify: pylint pep8

isort:
	isort --recursive pyvmd
