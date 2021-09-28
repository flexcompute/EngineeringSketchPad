#
test:
	python --version
	$(VALGRIND_COMMAND) python -u -m unittest discover

coverage:
	coverage run -m unittest discover
	coverage html
	open htmlcov/intex.html

clean:
	rm -rf __pycache__
	rm -f *.pyc

cleanall: clean
	rm -rf *UnitTest*
