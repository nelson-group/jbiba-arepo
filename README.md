# Arepo with Forced Isotropic Turbulence Module

This repository contains the public version of the Arepo code, along with additional code developed to run simulations of forced isotropic turbulence.

The primary extension can be found in:

    src/turbulence

## Compilation

When compiling the code, the following configuration options must be enabled:

    EXTERNALGRAVITY
    TURBULENCE

## Method

The forcing method implemented here follows the approach described in:

- Federrath et al. (2010), *Comparing the statistics of interstellar turbulence in simulations and observations: Solenoidal versus compressive turbulence forcing*. https://arxiv.org/abs/0905.1060  

An additional useful reference implementation by the original authors is available at:

- https://github.com/chfeder/turbulence_generator

---

The Arepo public version Readme:

Arepo public version
====================

AREPO is a massively parallel code for gravitational n-body 
systems and hydrodynamics, both on Newtonian as well as 
cosmological background. It is a flexible code that can be 
applied to a variety of different types of simulations, offering 
a number of sophisticated simulation algorithms. An description 
of the numerical algorithms employed by the code is given in the 
original code papers (Springel 2010, MNRAS, 401, 791; 
Pakmor et al. 2011, MNRAS, 418, 1392; Pakmor and Springel 2013, 
MNRAS, 432, 176; Pakmor et al. 2016, MNRAS,455,1134) and the 
release paper of this version (Weinberger et al. 2019). 

A user guide can be found under `/documentation`, which also 
includes a 'getting started' section, which is recommended for 
new users. An html version of the user guide can be created using
sphinx (https://www.sphinx-doc.org) by typing

    cd ./documentation/
    make html
    
and displayed by opening `./documentation/build/html/index.html`.

A full version of the user guide is also available on the Arepo 
homepage.
