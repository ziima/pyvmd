VMD = vmd

# All targets are phony
.PHONY: test coverage pylint flake8 pepify isort check-isort

test:
	${VMD} -python -dispdev none -e pyvmd/tests/__init__.py -args discover

coverage:
	python-coverage erase
	-rm -r htmlcov
	-COVERAGE=1 ${VMD} -python -dispdev none -e pyvmd/tests/__init__.py -args discover
	-python-coverage html -d htmlcov

pylint:
	-PYTHONPATH="/usr/lib/vmd/scripts/python" pylint pyvmd --reports=no

flake8:
	-flake8 pyvmd

pepify: pylint flake8

isort:
	isort --recursive pyvmd

check-isort:
	isort --check-only --diff --recursive pyvmd
