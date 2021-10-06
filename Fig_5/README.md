## Running simulations in Fig. 5
The simulation can be started from MATLAB by running the script
`matlab_Fig_5_run.m`.

If you are using macOS or are using Octave instead of
MATLAB, you will have to change the command to compile the `integrate_eqns.c`
at the top of the file.

By default, only 100 consolidation cycles are simulated to reduce simulation
time. In the manuscript, 2000 consolidation cycles were simulated, of which
the first 950 were used to initialize the weight matrices.

For plotting results from a previously run simulation, run the script
`matlab_Fig_5_plot.m`
