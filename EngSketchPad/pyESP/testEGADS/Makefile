#
test:
	python --version
	$(VALGRIND_COMMAND) python -u -m unittest discover

clean:
	rm -rf __pycache__
	rm -f *.pyc
