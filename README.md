# VAST Post-Processing

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

This repository contains the code of VAST Post-Processing, which further
processes ASKAP-VAST observation data.

## Features

* Code base in `Python 3.9.x`
* Crops and corrects FITS images
* Crops and corrects catalogues
* Generates MOCs and STMOCs

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

## Usage
This project is separated by module. The following modules are available by
command. 
- `link_neighbours`
- `convolve_neighbours`
- `swarp_neighbours`
- `selavy_prepare`
- `crop_fields`

To run the main program, run
```
vast_post_processing --data-root <DIRECTORY>
```

## Testing
The test suite expects at least one image observation, and its corresponding
data files. To check for the existence and organisation of the test data, run
the data test subsuite. 

1. Follow the installation instructions above, and if not already entered, enter
   the `vast-post-processing` environment with 
```
poetry shell
```
2. Test for local test data existence with 
```
pytest tests/data/test_data.py
```
3. If all tests pass, proceed to develop your own tests. Otherwise, your test
data directory is not initialised or configured correctly. Run the
`pull_data` script with 
```
pull_test_data
```
This process will take several minutes. 

Note this script requires `ssh` access to the `vast-data` virtual system. [Read
here for more
information](https://github.com/askap-vast/vast-project/wiki/Nimbus:-SSH-&-Downloading-Data). 

You can edit various settings and which observations to pull in
`tests/data/pull_config.yaml`.

4. After test data has been downloaded, the fixtures found in the various
`conftest.py` files will resolve, and all current tests should pass. Future
tests can now be implemented as well.

## Contributors

|Contributor|Organisation|Department|
|:---|:---|:---|
|[Dougal Dobie](https://github.com/ddobie)|[Swinburne University of Technology](https://www.swinburne.edu.au/research/our-research/access-our-research/find-a-researcher-or-supervisor/researcher-profile/?id=ddobie)|[Centre for Astrophysics and Supercomputing](https://www.swinburne.edu.au/research/centres-groups-clinics/centre-for-astrophysics-supercomputing/)|
|[Andrew O'Brien](https://github.com/marxide)|[University of Wisconsin-Milwaukee](https://uwm.edu/physics/people/obrien-andrew/)|[Center for Gravitation, Cosmology & Astrophysics](https://cgca.uwm.edu/people.html)|
|[Akash Anumarlapudi](https://github.com/AkashA98)|[University of Wisconsin-Milwaukee](https://uwm.edu/physics/people/anumarlapudi-akash/)|[Center for Gravitation, Cosmology & Astrophysics](https://cgca.uwm.edu/people.html)|
|[Mubdi Rahman](https://github.com/mubdi)|[Sidrat Research](https://sidratresearch.com/whoarewe.html)||
|[Hansen Jiang](https://github.com/hansenjiang)|[Sidrat Research](https://sidratresearch.com/whoarewe.html)||


<!-- ## Acknowledgements

The VAST Post-Processing development was supported by:

* Someone -->
