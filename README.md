# VAST Post-Processing

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

This repository contains the code of VAST Post-Processing, which further
processes ASKAP-VAST observation data.

## Features

* Code base in `Python 3.9.x`
* Crops and corrects images
* Generates MOCs

## Installation
1. Clone this repository onto your local machine. 
```
git clone git@github.com:askap-vast/vast-post-processing.git
```
2. Navigate to the repository. 
```
cd vast-post-processing/
```
3. Install [Poetry](https://python-poetry.org/) and [Poe the
   Poet](https://poethepoet.natn.io/).
```
pip install poetry
pip install poethepoet
```
4. Install the rest of this package's dependencies with Poetry.
```
poetry install
```
5. The package is now installed. You can enter a virtual environment with 
```
poetry shell
```

## Testing
# TODO WRITE
1. Download data

<!-- ## Usage
1. To be written -->

## Contributors

* Andrew O'Brien – [Department of Physics, University of Wisconsin-Milwaukee](https://uwm.edu/physics/research/astronomy-gravitation-cosmology/)
* Dougal Dobie – [Centre for Astrophysics and Supercomputing, Swinburne
  University of
  Technology](https://www.swinburne.edu.au/research/our-research/access-our-research/find-a-researcher-or-supervisor/researcher-profile/?id=ddobie)
* Akash Anumarlapudi - [Department of Physics, University of
  Wisconsin-Milwaukee](https://uwm.edu/physics/people/anumarlapudi-akash/) 
* Mubdi Rahman - [Sidrat Research](https://www.sidratresearch.com/)
* Hansen Jiang - [Sidrat Research](https://www.sidratresearch.com/)

<!-- ## Acknowledgements

The VAST Post-Processing development was supported by:

* Someone -->
