SHELL = /bin/bash -e 

MY_NOSE_FLAGS?=-v -s
MY_CRAM_FLAGS?=-v
THISDIR:=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))


test-fast: pylint

wheel:
	# This is basically a syntax check.
	python setup.py bdist_wheel
build:
	pip install -v --user --edit .
	rm -f *.xml
pylint:
	pylint --errors-only pbtranscript_internal_validation/
install:
	python setup.py install
test:
	# Unit tests
	#find tests/unit -name "*.py" |grep -v test_validate_nti.py | xargs nosetests
	nosetests -s tests/unit/test_functions.py
#slowtest:
#	nosetests -s tests/unit/test_validate_nti.py
clean:
	rm -rf dist/ build/ *.egg-info
	rm -rf doc/_build
	find . -name "*.pyc" | xargs rm -f
	rm -rf dist/ LOCAL/
pip-install:
	@which pip > /dev/null
	@pip freeze|grep 'pbtranscript-internal-validation=='>/dev/null \
      && ( pip uninstall -y pbtranscript-internal-validation) \
      || true
	@python setup.py build
	@pip install --no-index \
          --install-option="--install-scripts=$(PREFIX)/bin" \
          ./
push:
	git push origin HEAD:feature/SE-1072-unify-isoseq-isoseq2-validation

.PHONY: all build install test-fast clean pip-install
