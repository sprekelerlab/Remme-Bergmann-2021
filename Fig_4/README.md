## Running simulations in Fig. 4
The simulation can be started from MATLAB by running the script
`matlab_Fig_4_run.m`.

If you are using macOS, you might have to change the command to compile the
`integrate_eqns.c` at the top of `matlab_Fig_4_run.m`

Depending on whether you use MATLAB or Octave, you need to change the seed
function in line 61/62 of `matlab_Fig_4_run.m` (the current default assumes
Octave is being used).

By default, no lesion is simulated. To simulate a lesion, change the
`PP_lesion` parameter in `matlab_Fig_4_run.m`.

For plotting results from a previously run simulation, run the script
`matlab_Fig_4_plot.m`. Note, that the `N_cycle` parameter in the
`matlab_Fig_4_plot.m` needs to be the same as in `matlab_Fig_4_run.m`.
