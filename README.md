# MCRaT
fully parallelize MCRaT
changed mc.par to include frame injections for small and large angles, as well as independent weights and radii of injection.

Also included min and max nmber of photons so the program can automatically adjust the photon weight to produce photons within min and max. This leads to better statistics.

Can now specify number of threads. Half goes to splitting up simulation in angle and other half goes to splitting up angles via injection frames

Added log files for each thread to see progress
