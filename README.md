# boxfit

The BOXFIT gamma-ray burst afterglow fit and light curve generator code is a numerical implementation of the work described in [Van Eerten+ 2012](https://ui.adsabs.harvard.edu/abs/2012ApJ...749...44V/abstract). The code is capable of calculating light curves and spectra for arbitrary observer times and frequencies and of performing (broadband) data fits using the downhill simplex method combined with simulated annealing. The flux value for a given observer time and frequency is a function of various variables that set the explosion physics (energy of the explosion, circumburst number density and jet collimation angle), the radiative process (magnetic field generation efficiency, electron shock-acceleration efficiency and synchrotron power slope for the electron energy distribution) and observer position (distance, redshift and angle).

The dynamics of the afterglow blast wave have been calculated in a series of 114 high-resolution two-dimensional jet simulations performed with the RAM adaptive-mesh refinement relativistic hydrodynamics (RHD) code [Zhang+ 2006](https://ui.adsabs.harvard.edu/abs/2006ApJS..164..255Z/abstract). The results of these calculations have been compressed and stored in a series of `box' data files and BOXFIT calculates the fluid state for arbitrary fluid variables using interpolations between the data files and analytical scaling relations. End-users of BOXFIT do not need to perform RHD simulations themselves.

The code can be run both in parallel and on a single core. Because a data fit takes many iterations, this is best done in parallel. Single light curves and spectra can readily be done on a single core. Use of the code is completely free, but we request that [Van Eerten+ 2012](https://ui.adsabs.harvard.edu/abs/2012ApJ...749...44V/abstract) be cited in scientific publications using the fitting algorithms from BOXFIT. This also applies to derived code, plus possible additional citations depending on the modifications.

## Getting Started

A more detailed guide to working with the code is provided by the guide (pdf) included with installation.

### Prerequisites

BOXFIT makes use of the [hdf5](https://www.hdfgroup.org/solutions/hdf5/) file format, and requires the developer libraries for hdf5 to be installed on the system. Parallel BOXFIT requires MPI libraries to be installed.

### Installing

First clone the repository:
```
git clone https://github.com/hveerten/boxfit
```
In the */src* subdirectory, set up a makefile starting from the *makefile.template* file:
```
cd boxfit/src
cp makefile.template makefile
```


### Running the code

## Authors

* **Hendrik van Eerten** - *University of Bath* - [github](https://github.com/hveerten) - [University of Bath](https://researchportal.bath.ac.uk/en/persons/hendrik-van-eerten)

