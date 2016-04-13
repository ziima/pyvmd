VMD = vmd

# All targets are phony
.PHONY: test coverage pylint pepify isort check-isort check-flake8

test:
	${VMD} -python -dispdev none -e pyvmd/tests/__init__.py -args discover

coverage:
	python-coverage erase
	-rm -r htmlcov
	-COVERAGE=1 ${VMD} -python -dispdev none -e pyvmd/tests/__init__.py -args discover
	-python-coverage html -d htmlcov

pylint:
	-PYTHONPATH="/usr/lib/vmd/scripts/python" pylint pyvmd --reports=no

pepify: pylint check-flake8

isort:
	isort --recursive pyvmd

check-isort:
	isort --check-only --diff --recursive pyvmd

check-flake8:
	flake8 --format=pylint pyvmd
