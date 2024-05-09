#
test:
	python --version
	$(VALGRIND_COMMAND) python -u -m unittest discover --verbose

clean:
	rm -rf __pycache__
	rm -f *.pyc

cleanall: clean
	rm -rf workDir_* tmp_OpenCSM*
