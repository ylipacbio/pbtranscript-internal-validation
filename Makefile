SHELL = /bin/bash -e 

all: check build install

check:
	pylint --errors-only pbtranscript_internal_validation/*.py

build:
	python setup.py build --executable="/usr/bin/env python"

bdist:
	python setup.py build --executable="/usr/bin/env python"
	python setup.py bdist --formats=egg
pylint:
	pylint --errors-only pbtranscript_internal_validation/
install:
	python setup.py install

develop:
	python setup.py develop

test:
	# Unit tests
	find tests/unit -name "*.py" |grep -v test_validate_nti.py | xargs nosetests
	# End-to-end tests
	#find tests/cram -name "*.t" | xargs cram

slowtest:
	nosetests -s tests/unit/test_validate_nti.py

doc:
	sphinx-apidoc -T -f -o doc src/ && cd doc && make html
docs: doc

clean: doc-clean
	rm -rf dist/ build/ *.egg-info
	rm -rf doc/_build
	find . -name "*.pyc" | xargs rm -f
	rm -rf dist/

doc-clean:
	rm -f doc/*.html

pip-install:
	@which pip > /dev/null
	@pip freeze|grep 'pbtranscript-internal-validation=='>/dev/null \
      && ( pip uninstall -y pbtranscript-internal-validation) \
      || true
	@python setup.py build
	@pip install --no-index \
          --install-option="--install-scripts=$(PREFIX)/bin" \
          ./

.PHONY: all build bdist install develop test doc clean pip-install
