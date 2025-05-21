# Documentation

This directory contains the documentation for this project. 

The modules of this project are documented using docstrings
in [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html) style and
comments, and formatted using code blocking, parentheses, and other
[PEP8](https://peps.python.org/pep-0008/) style guidelines, using [the Black
formatter](https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html). 

They are auto-generated in `HTML` format with [`sphinx`](https://www.sphinx-doc.org/en/master/index.html). 

## Included

1. `source/`
    1. `conf.py`
    2. `index.rst`
    3. `modules.rst`
    4. `vast_post_processing.rst`
2. `make.bat`
3. `Makefile`

## Instructions

To view `sphinx` documentation for this project, navigate to the root of the
package (i.e. `vast-post-processing`), and
run
```
poetry install
poetry shell
cd docs
make html
```
The pages are built in `vast-post-processing/docs/build/html`, and you can load
the index page by opening `index.html` in a browser. 
