## Running simulations in Fig. 3
To run this code, the [NEURON](https://neuron.yale.edu/neuron/) simulator
needs to be installed.

### Prepare the NEURON simulation
Before running the code, the `.mod` files in the directory `_mod` need to be
compiled with the application `nrnivmodl` (or with `mknrndll`) that is included
with NEURON.
This will create the executable `_mod/<arch>/special`,
where `<arch>` depends on your CPU architecture, e.g. `x86_64` for 64 bit
systems. This executable should be copied manually to the same directory as
`neuron_Fig_3.hoc`. Here are the commands to be executed from inside the
`Fig_3` directory:
```bash
cd _mod
nrnivmodl *.mod
# Check what the <arch> directory is called on your system, e.g. x86_64 and
# modify the next command accordingly
cp <arch>/special ..
```
Now you should have the `special` binary located in the same folder as
`neuron_Fig_3.hoc`.


# Run the simulation
The simulation can be started from MATLAB by running the script
`matlab_Fig_3_run.m`. If you are using macOS, you will have to change the
command to run NEURON in line 95.

For plotting results from a previously run simulation, run the script
`matlab_Fig_3_plot.m`
