init-linux:
	conda env create --file=linux-environment.yml

init-windows:
	conda env create --file=windows-environment.yml

compile:
	python setup.py

requirements-linux:
	conda env export > linux-environment.yml --no-builds

requirements-windows:
	conda env export > windows-environment.yml --no-builds

reformat:
	autoflake --in-place --remove-all-unused-imports  --ignore-init-module-imports  --remove-duplicate-keys \
	--remove-unused-variables -r .
	isort .
	yapf --style "google" -vv -ir .

test:
	python validate/validate.py
