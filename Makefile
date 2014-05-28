VMD = vmd

test:
	${VMD} -python -dispdev none -e pyvmd/tests/__init__.py -args discover

coverage:
	python-coverage erase
	-rm -r htmlcov
	COVERAGE=1 ${VMD} -python -dispdev none -e pyvmd/tests/__init__.py -args discover
	python-coverage html -d htmlcov

pylint:
	PYTHONPATH="/usr/lib/vmd/scripts/python" pylint pyvmd
