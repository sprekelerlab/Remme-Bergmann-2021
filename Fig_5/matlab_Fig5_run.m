% Memory consolidation - Figure 5
% Michiel Remme, May 2016
% memory consolidation in a hierarchical network
% for analyzing and plotting data after running this simulation, use matlab_Fig_5_plot.m

clear all

% the gnu scientific library (gsl) needs to be installed because of the random number generator used in the c-code

% compile the c-code
% on mac
% mex -I/usr/local/include -L/usr/local/lib -lgsl -lm integrate_eqns.c
% on linux machines need to add lgslcblas
mex -lgsl -lgslcblas -lm integrate_eqns.c

%% Parameters

seed_init   = 100;
N_cycle     = 100;  % (2000) number of consolidation cycles (article used 2000 cycles of which the first 950 were used for initialization of the weight matrices)
nsec        = 150;  % (150) sec

N_layer     = 8;    % (8) Number of input layers = number of output layers
N_cell      = 256;  % (256) No cells per layer

dt          = 5*1e-3; % sec
t_del       = 5*1e-3; % sec - transmission delay between layers
it_del      = 1+t_del/dt; % 1 + delay in # iterations
N_HPC_init_cycle = 30; % number of memories stored in HPC before running actual initialization simulations

A_mean      = 10;   % 10 (Hz) mean input layer activity
A_std       = 5;    % 5 (Hz) stdev input layer activity

alpha       = 1.00008; % ratio of LTD vs LTP
tau_LTP     = 20*1e-3; % sec
tau_LTD     = 20*1e-3; % sec

w_max       = 2/N_cell;
lambda      = 0.2*w_max; % learning rate
lambda_q    = 0.5;  % lambda decreases with factor lambda_q per layer
HPC_mem_fac = 0.5;  % relative contribution of new memories

%% %%%%%%%%%%%%%%%%%%%%%%%%%
% run batch

mkdir('_results')

seed    = seed_init;
#rng(seed); % in matlab
rand('state', seed); % in octave

%%% INITIALIZE HPC MATRIX
W_HPC_init = binornd(1,0.5,N_cell,N_cell);
W_HPC_init = bsxfun(@rdivide,W_HPC_init,sum(W_HPC_init,2)); % L1 normalize all rows of W_HPC_init
for i_cycle = 1:N_HPC_init_cycle
    W_HPC_mem       = binornd(1,0.5,N_cell,N_cell);
    W_HPC_mem       = bsxfun(@rdivide,W_HPC_mem,sum(W_HPC_mem,2)); % L1 normalize all rows of W_HPC_update
    W_HPC_init      = W_HPC_init + W_HPC_mem*HPC_mem_fac/(1-HPC_mem_fac);
    W_HPC_init      = bsxfun(@rdivide,W_HPC_init,sum(W_HPC_init,2)); % L1 normalize all rows of W_HPC_init
end

W_init = w_max*rand(N_cell,N_cell,N_layer);

file = sprintf('_results/data_Fig5_seed_%d_Nlayer_%d_cycle_0',seed_init,N_layer);
save(file); % save all parameters in this file - W_HPC_mem is the memory that will be tracked and that has not yet been learned

%%% RUN CONSOLIDATION CYCLES
disp 'running consolidation cycles...'
for i_cycle = 1:N_cycle
    disp(i_cycle)
    % run consolidation cycle
    params  = [nsec dt N_layer N_cell A_mean A_std it_del w_max lambda lambda_q alpha tau_LTP tau_LTD seed];
    [W_,A_]         = integrate_eqns(W_HPC_init,W_init,params);
    % save weight matrices
    file = sprintf('_results/data_Fig5_seed_%d_Nlayer_%d_cycle_%d',seed_init,N_layer,i_cycle);
    save(file,'W_HPC_mem','W_HPC_init','W_','A_');
    % weight matrices update
    W_HPC_mem       = binornd(1,0.5,N_cell,N_cell);
    W_HPC_mem       = bsxfun(@rdivide,W_HPC_mem,sum(W_HPC_mem,2)); % L1 normalize all rows of W_HPC_mem
    W_HPC_init      = W_HPC_init + W_HPC_mem*HPC_mem_fac/(1-HPC_mem_fac);
    W_HPC_init      = bsxfun(@rdivide,W_HPC_init,sum(W_HPC_init,2)); % L1 normalize all rows of W_HPC_init
    W_init          = W_;
    seed            = seed + 1;
end

%%
