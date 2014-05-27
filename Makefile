VMD = vmd

test:
	${VMD} -python -dispdev none -e pyvmd/tests/__init__.py -args "discover"

pylint:
	PYTHONPATH="/usr/lib/vmd/scripts/python" pylint pyvmd
