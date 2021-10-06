## Running simulations in Fig. 3
To run this code, the NEURON simulator needs to be installed.

Before running the code, the MOD files in the directory `_mod` need to be
compiled with the application `mknrndll` (or with `nrnivmodl`) that is included
with NEURON. This will create the executable `special` in a subdirectory of
_mod. This executable should be copied manually to the same directory as
`neuron_Fig_3.hoc`.

The simulation can be started from MATLAB by running the script
`matlab_Fig_3_run.m`. If you are using macOS, you will have to change the
command to run NEURON in line 95.

For plotting results from a previously run simulation, run the script
`matlab_Fig_3_plot.m`
