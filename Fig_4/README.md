## Running simulations in Fig. 4
The simulation can be started from MATLAB by running the script
`matlab_Fig_4_run.m`.

If you are using macOS or are using Octave instead of
MATLAB, you will have to change the command to compile the `integrate_eqns.c`
at the top of the file.

By default, no lesion is simulated. To simulate a lesion, change the
`PP_lesion` parameter in `matlab_Fig_4_run.m`.

For plotting results from a previously run simulation, run the script
`matlab_Fig_4_plot.m`
