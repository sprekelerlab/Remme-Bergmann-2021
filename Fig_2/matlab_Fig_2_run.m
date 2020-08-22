% Memory consolidation - Figure 2
% Michiel Remme, May 2016
% run CA1 pyramidal cell model in NEURON from Matlab and plot the results
% note, for only plotting results from previously run simulation, use 

clear all

%% 1) set simulution parameters and save to NEURON-readable file

tstop       = 5*60*1000; % ms; Note, should not exceed 'tend' used in create_spike_times.m
dt          = 0.025; % ms
input_vec_dt= 0.05; % ms; Note, should match 'dt' used in create_spike_times.m
rec_start   = 10000; % ms - record voltage from this time...
rec_stop    = 11000; % ms - until this time...
rec_dt      = 0.025; % ms - with this time step
save_we_ivl = 1; % sec - interval for saving weight matrices
T_cycle     = 0.5; % sec

Ne_sc       = 2500; % number of SC inputs
Ne_pp       = 500; % number of PP inputs

vth_sc      = -30; % (mV) threshold for stdp learning rule and for recording spike times
vth_pp      = -30; % (mV) 

we_sc_max   = 0.0004; % uS
we_pp_max 	= 0.00014; % uS
we_sc_min	= 0;
we_pp_min	= 0;

del_sc      = 5; % ms
del_pp      = 0; % ms
lambda_sc   = 0*0.05; % relative weight changes
lambda_pp   = 1*0.05; % relative weight changes
alpha       = 1.05; % fraction depression versus potentiation

dir_init = '_initialization_files';
dir_lib  = '_lib';
dir_res  = '_results/';

% Generate parameters file for Neuron
fid = fopen(fullfile(dir_init,'params.m'), 'wt');
fprintf(fid,'tstop = %g\n',tstop);
fprintf(fid,'dt = %g\n',dt);
fprintf(fid,'rec_dt = %g\n',rec_dt);
fprintf(fid,'rec_start = %g\n',rec_start);
fprintf(fid,'rec_stop = %g\n',rec_stop);
fprintf(fid,'Ne_sc = %g\n',Ne_sc);
fprintf(fid,'Ne_pp = %g\n',Ne_pp);
fprintf(fid,'vth_sc = %g\n',vth_sc);
fprintf(fid,'vth_pp = %g\n',vth_pp);
fprintf(fid,'we_sc_max = %g\n',we_sc_max);
fprintf(fid,'we_sc_min = %g\n',we_sc_min);
fprintf(fid,'we_pp_max = %g\n',we_pp_max);
fprintf(fid,'we_pp_min = %g\n',we_pp_min);
fprintf(fid,'lambda_sc = %g\n',lambda_sc);
fprintf(fid,'lambda_pp = %g\n',lambda_pp);
fprintf(fid,'alpha = %g\n',alpha);
fprintf(fid,'del_sc = %g\n',del_sc);
fprintf(fid,'del_pp = %g\n',del_pp);
fprintf(fid,'save_we_ivl = %g\n',save_we_ivl);
fclose(fid);

%% 2) generate spike times and place/grid fields and save to files for NEURON

run(fullfile('_lib','create_spike_times.m'))

%% 3) load file which fibers to potentiate and save initial weights to Neuron-readable files

% create initial weight vectors
% phase of ca1 'place field'
ca1_phase = 0.5*T_cycle*1000/input_vec_dt; % time in cycle

% sc weight vectors (initialize with sc input generating a single place field in the CA1 cell)
sc_pot_vec = (sc_phase<ca1_phase) & (sc_phase + T_cycle*1000/input_vec_dt./(sc_spacing_vec*2))>ca1_phase;
we_sc_init  = we_sc_min*ones(Ne_sc,1);
we_sc_ind = sc_pot_vec;
we_sc_init(we_sc_ind) = we_sc_max;

% pp weight vectors
we_pp_init  = we_pp_min*ones(Ne_pp,1);
we_pp_ind   = randperm(Ne_pp,0.5*Ne_pp); % potentiate random 50% of pp fibers
we_pp_init(we_pp_ind) = we_pp_max;

% save initial weights
fid = fopen(fullfile(dir_init,'we_sc_init.dat'), 'wt');
fprintf(fid,'%g\n',we_sc_init);
fclose(fid);
fid = fopen(fullfile(dir_init,'we_pp_init.dat'), 'wt');
fprintf(fid,'%g\n',we_pp_init);
fclose(fid);

%% 4) run Neuron

system('./special -nobanner -nogui neuron_Fig_2.hoc'); % for mac
% system('LD_LIBRARY_PATH="" ./special -nobanner -nogui neuron_ca1_consolidation.hoc'); % for LINUX

%% 5) plot data

run(fullfile(dir_lib,'plot_results'))
