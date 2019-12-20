# boxfit

The BOXFIT gamma-ray burst afterglow fit and light curve generator code is a numerical implementation of the work described in Van Eerten+ 2012. The code is capable of calculating light curves and spectra for arbitrary observer times and frequencies and of performing (broadband) data fits using the downhill simplex method combined with simulated annealing. The flux value for a given observer time and frequency is a function of various variables that set the explosion physics (energy of the explosion, circumburst number density and jet collimation angle), the radiative process (magnetic field generation efficiency, electron shock-acceleration efficiency and synchrotron power slope for the electron energy distribution) and observer position (distance, redshift and angle).

The dynamics of the afterglow blast wave have been calculated in a series of 114 high-resolution two-dimensional jet simulations performed with the RAM adaptive-mesh refinement relativistic hydrodynamics (RHD) code (Zhang+ 2006). The results of these calculations have been compressed and stored in a series of `box' data files and BOXFIT calculates the fluid state for arbitrary fluid variables using interpolations between the data files and analytical scaling relations. End-users of BOXFIT do not need to perform RHD simulations themselves.

The code can be run both in parallel and on a single core. Because a data fit takes many iterations, this is best done in parallel. Single light curves and spectra can readily be done on a single core. Use of the code is completely free, but we request that Van Eerten+ 2012 be cited in scientific publications using the fitting algorithms from BOXFIT. This also applies to derived code, plus possible additional citations depending on the modifications.

## Getting Started

## Authors

* **Hendrik van Eerten** - *Initial work* - [test](https://github.com/hveerten)

