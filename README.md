# MCRaT

Added tests for Serial Cases and for parallel cases, also found bug with lorentx boost where i was calculating beta incorrectly and other bug where I was not passing singleElectron function the 4 momentum of the photon in the comoving frame.

Found another bug where calling rand for muliple threads calling it gives undefined behavior.

Solved bug for multiple threads calling rand by creating an array of random number generators for each thread.

