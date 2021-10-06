## Running simulations in Fig. 5
The simulation can be started from MATLAB by running the script
`matlab_Fig_5_run.m`.

If you are using macOS, you might have to change the command to compile the
`integrate_eqns.c` at the top of `matlab_Fig_4_run.m`

Depending on whether you use MATLAB or Octave, you need to change the seed
function in line 61/62 of `matlab_Fig_4_run.m` (the current default assumes
Octave is being used).

By default, only 100 consolidation cycles are simulated to reduce simulation
time. In the manuscript, 2000 consolidation cycles were simulated, of which
the first 950 were used to initialize the weight matrices. You can change this
via the `N_cycle` parameter.

For plotting results from a previously run simulation, run the script
`matlab_Fig_5_plot.m`
