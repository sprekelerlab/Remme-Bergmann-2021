## Running simulations in Fig. 1
The simulation can be started from Octave or MATLAB by running the script
`matlab_Fig_1.m`.

If you are using macOS instead of Linux, you might have to change the command
to compile `integrate_eqns.c` at the top of the `matlab_Fig_1.m` file.

To turn plasticity off in the PP and on in the SC (Fig. 1F), set `pp_plastic = false;`
in `matlab_Fig_1.m`.
