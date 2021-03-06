This repository provides the code to reproduce the simulations in the paper

[**Hebbian plasticity in parallel synaptic pathways: A circuit mechanism for systems memory consolidation**](https://doi.org/10.1371/journal.pcbi.1009681)

*Michiel Remme\*, Urs Bergmann\*, Denis Alevi, Susanne Schreiber, Henning Sprekeler+, Richard Kempter+*
____

To run the code in this repository, you need to install the following software:
- [GNU Octave](https://www.gnu.org/software/octave/index) or [MATLAB](https://de.mathworks.com/products/matlab.html)
- The [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) (for
  the random number generation in the simulations of Fig. 1, 4 and 5)
- The [NEURON](https://neuron.yale.edu/neuron/) simulator (for the simulations in Fig. 3.)

Detailed instructions on how to run the simulations can be found in the README files of the `Fig` folders.

This code was tested on macOS and Linux. Depending on whether MATLAB or Octave
is used and whether macOS or Linux is used, some files have to be modified as
described in the separate `Fig` folders. All code was last run with Octave
6.3.0, GSL 2.7 and NEURON 7.8.1.
