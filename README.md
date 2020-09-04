# Welcome to the Github page for the Monte Carlo Radiation Transfer Code 

The Monte Carlo Radiation Transfer (MCRaT; pronounced \textit{Em-Cee-Rat}) code is a next generation radiation transfer code that can be used to analyze the radiation signature expected from astrophysical outflows. The code is written in C and uses the Message Passing Interface (MPI) library for inter-process communication, the Open Multi-Processor (OpenMP) library for intra-node communication, the GNU Scientific Library, and the HDF5 library for parallel I/O.

MCRaT injects photons in a FLASH simulation and individually propagates and compton scatters the photons through the fluid until the end of the simulation. This process of injection and propagating occurs for a user specified number of times until there are no more photons to be injected. Users can then construct light curves and spectra with the MCRaT calculated results. The hydrodynamic simulations used with this version of MCRaT must be in 2D, however, the photon propagation and scattering is done in 3D by assuming cylindrical symmetry. Additionally, MCRaT uses the full Klein–Nishina cross section including the effects of polarization, which can be fully simulated in the code.

Currently, MCRaT works with FLASH hydrodynamic simulations and PLUTO AMR simulations, with both 2D spherical (r, ![equation](https://latex.codecogs.com/gif.latex?%5Ctheta)) and 2D cartesian ((x,y) and (r,z)).
<!-- for tex: https://stackoverflow.com/questions/35498525/latex-rendering-in-readme-md-on-github -->

MCRaT is also able to emit and absorb thermal synchrotron photons and keep track of them as the photons, and the outflow, propagates. 

The video below shows an example of MCRaT scattering photons through a Gamma Ray Burst jet. MCRaT allows us to track the evolution fo the photon spectrum as they photons propagate through the jet and as the jet propagates through space. Additionally, we can keep track of how well the photons are in equilibrium with the matter in the jet, as indicated by the effective temperatures. 

[![](https://img.youtube.com/vi/pjkAyGUsJro/0.jpg)](https://www.youtube.com/watch?v=pjkAyGUsJro)
 

The code was initially written in python by Dr. Davide Lazzati as a proof of concept. The code was then translated into C by Tyler Parsotan and made to use the OpenMP, MPI, HDF5, and GNU Scientific libraries. MCRaT is highly parallelized and is easy to set up and use.

There are also python files with documented functions to process the MCRaT output located at: https://github.com/parsotat/ProcessMCRaT.

In using MCRaT and the ProcessMCRaT codes, we ask that you cite the following papers: [Lazzati (2016)](https://doi.org/10.3847/0004-637X/829/2/76); [Parsotan & Lazzati (2018)](https://doi.org/10.3847/1538-4357/aaa087); [Parsotan et al. (2018)](https://doi.org/10.3847/1538-4357/aaeed1); [Parsotan et. al. (2020)](https://doi.org/10.3847/1538-4357/ab910f).
