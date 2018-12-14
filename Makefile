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
test:
	# Unit tests
	#find tests/unit -name "*.py" |grep -v test_validate_nti.py | xargs nosetests
	nosetests -s tests/unit/test_functions.py
slowtest:
	nosetests -s tests/unit/test_validate_nti.py
pip-install:
	@which pip > /dev/null
	@pip freeze|grep 'pbtranscript-internal-validation=='>/dev/null \
      && ( pip uninstall -y pbtranscript-internal-validation) \
      || true
	@python setup.py build
	@pip install --no-index \
          --install-option="--install-scripts=$(PREFIX)/bin" \
          ./
autopep8:
	find pbtranscript_internal_validation/ -type f -name '*.py' | xargs autopep8  --in-place --max-line-length 120 -ir -j 0
	find setup.py -type f -name '*.py' | xargs autopep8  --in-place --max-line-length 120 -ir -j 0
autofmt:
	find pbtranscript_internal_validation/ -type f -name '*.py' | xargs autoflake --in-place --remove-unused-variables --expand-star-imports
	find pbtranscript_internal_validation/ -type f -name '*.py' | xargs autopep8  --in-place --max-line-length 120 -ir -j 0
	find setup.py -type f -name '*.py' | xargs autopep8  --in-place --max-line-length 120 -ir -j 0
vulture:
	vulture isocollapse/

clean:
	find . -type f \( -name '*.pyc' -o -name '*.so' -o -name '*.o' -o -name '*.cpp' \) -delete
	rm -rf dist/ build/ *.egg-info doc/_build LOCAL/

.PHONY: all build install test-fast clean pip-install
