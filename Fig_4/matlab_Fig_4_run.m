% Memory consolidation - Figure 4
% Michiel Remme, May 2016
% memory consolidation in a hippocampal network
% for analyzing and plotting results, run the script matlab_Fig_4_plot.m

clear all

% the gnu scientific library (gsl) needs to be installed because of the random number generator used in the c-code

% compile the c-code
% on mac
% mex -I/usr/local/include -L/usr/local/lib -lgsl -lm integrate_eqns.c
% on linux machines need to add lgslcblas
mex -lgsl -lgslcblas -lm integrate_eqns.c

%% Parameters
PP_lesion   = 0; % 0 = no lesion // 1 = lesion before training // 2 = lesion after 21 cycles
N_trial     = 1; % (10) // ~ an hour simulation time per trial
N_cycle     = 31; % (31) number of consolidation cycles
nsec        = 150; % (150) sec
dt          = 5*1e-3; % sec
seed_range  = 1:N_trial; % random seeds used for each trial

% Number of cells per position or object population (EC, CA3, CA1)
N           = 16^2; % choose square of integer
env_len     = 1; % (meter) side of square box
N_obj       = 128; % number of objects

alpha       = 1.00025; % ratio of LTD vs LTP
tau_LTP     = 20*1e-3; % sec
tau_LTD     = 20*1e-3; % sec
A_EC_max    = 10; % max activity level
A_CA3_max   = 10;

W_PP_max    = 2/(2*N); % max weight of PP fibers
W_SC_id     = 1/4;
W_PP_init_scale = 1; %
lambda_PP   = 0.05*W_PP_max; % learning rate
t_del_CA3   = 5*1e-3; % sec - delay from EC to CA1 via SC
it_del_CA3  = 1+t_del_CA3/dt; % 1 + delay in # iterations (size of history array)

W_PPS_max   = W_PP_max; % PP to subiculum
W_SCS_id    = 1/2;
W_PPS_init_scale = 1; %
lambda_PPS  = 0.5*lambda_PP; %
t_del_CA1   = 1*dt; % sec - delay from CA1 to SUB
it_del_CA1  = 1+t_del_CA1/dt; % 1 + delay in # iterations

SC_update_fac   = 0.6; % overwriting each cycle in SC pathway
W_SC_po_init    = zeros(N,N);
N_mem           = 125;
N_init_cycles   = 0; % number of consolidation cycles to initialize WPP matrices

dir_res  = '_results/';
mkdir(dir_res)

%% run batch
for i_trial = 1:N_trial
    i_trial
    seed = seed_range(i_trial);
    #rng(seed); % in matlab
    rand('state', seed); % in octave

    %%% SET CA3 PLACE FIELD PARAMETERS AND OBJECT CODING ACTIVITIES
    pos_vec     = -env_len/2+env_len*(0.5:1:sqrt(N))./sqrt(N);
    [xvec yvec] = meshgrid(pos_vec);
    CA3_p_pf    = [xvec(:) yvec(:)];
    CA3_pf_var  = 0.1^2; % (cm)
    CA3_o_pf    = CA3_p_pf;
    CA3_o_A     = zeros(N,N_obj);
    for i=1:N_obj
        o_x   = -env_len/2+env_len*rand;
        o_y   = -env_len/2+env_len*rand;
        A_tmp = A_CA3_max*exp(-((o_x-CA3_o_pf(:,1)).^2 + (o_y-CA3_o_pf(:,2)).^2)./(2*CA3_pf_var));
        CA3_o_A(:,i) = A_tmp(randperm(N));
    end

    %%% SET EC GRID FIELD PARAMETERS AND OBJECT CODING ACTIVITIES
    gf_x    = env_len*rand(N,1);
    gf_y    = env_len*rand(N,1);
    gf_ori  = rand(N,1)*2/3*pi; % gf orientation (between 0-120 deg)
    gf_per  = linspace(2,6,N)'*2*pi/(env_len); % between 2 and 6 fields in environment
    EC_p_gf = [gf_x gf_y gf_ori gf_per];
    gf_x    = env_len*rand(N,1);
    gf_y    = env_len*rand(N,1);
    gf_ori  = rand(N,1)*2/3*pi; % gf orientation (between 0-120 deg)
    gf_per  = linspace(2,6,N)'*2*pi/(env_len); % between 2 and 6 fields in environment
    EC_o_gf = [gf_x gf_y gf_ori gf_per];
    EC_o_A  = zeros(N,N_obj);
    for i=1:N_obj
        o_x   = -env_len/2+env_len*rand;
        o_y   = -env_len/2+env_len*rand;
        A_tmp = 0;
        A_tmp =         cos((cos(        EC_o_gf(:,3)).*(o_x-EC_o_gf(:,1)) + sin(       EC_o_gf(:,3)).*(o_y-EC_o_gf(:,2))).*EC_o_gf(:,4));
        A_tmp = A_tmp + cos((cos(pi/3+   EC_o_gf(:,3)).*(o_x-EC_o_gf(:,1)) + sin(pi/3+  EC_o_gf(:,3)).*(o_y-EC_o_gf(:,2))).*EC_o_gf(:,4));
        A_tmp = A_tmp + cos((cos(2*pi/3+ EC_o_gf(:,3)).*(o_x-EC_o_gf(:,1)) + sin(2*pi/3+EC_o_gf(:,3)).*(o_y-EC_o_gf(:,2))).*EC_o_gf(:,4));
        A_tmp = A_EC_max*(A_tmp+1.5)/4.5;
        EC_o_A(:,i) = A_tmp(randperm(N));
    end

    %%% INITIALIZE ALL WEIGHT MATRICES
    %%% INITIALIZE SC WEIGHTS
    % update weights of sc cross fibers
    W_SC_po_init = zeros(N,N);
    for i_mem = 1:N_mem
        % location of object in circular environment
        p_ang  = rand*2*pi;
        p_dis  = rand + rand;
        if p_dis>1, p_dis = 2-p_dis; end
        p_x  = env_len/2*p_dis*cos(p_ang);
        p_y  = env_len/2*p_dis*sin(p_ang);
        % object identity
        o_ind  = randi(N_obj); % random object identity
        % SC update
        A_CA1_p_tmp     = A_CA3_max*exp(-((p_x-CA3_p_pf(:,1)).^2 + (p_y-CA3_p_pf(:,2)).^2)./(2*CA3_pf_var));
        W_SC_po_update  = A_CA1_p_tmp(randperm(N)) * CA3_o_A(:,o_ind)'; % REMAP CA1_p PATTERN
        W_SC_po_update  = W_SC_po_update./max(sum(W_SC_po_update,2));
        W_SC_po_init    = W_SC_po_init + W_SC_po_update*SC_update_fac/(1-SC_update_fac);
        W_SC_po_init    = W_SC_po_init./max(sum(W_SC_po_init,2));
    end
    %%% INITIALIZE PP + PPS WEIGHTS
%   W_PP_po_init    = W_PP_max*W_PP_init_scale*rand(N,N);
%   W_PP_pp_init    = W_PP_max*W_PP_init_scale*rand(N,N);
%   W_PPS_po_init   = W_PPS_max*W_PPS_init_scale*rand(N,N);
%   W_PPS_pp_init   = W_PPS_max*W_PPS_init_scale*rand(N,N);
    load W_PP_W_PPS_init % load weight matrices to avoid initialization artefacts (otherwise use previous 4 lines instead)
    W_PP_po_init    = W_PP_po_end(randperm(N*N));
    W_PP_pp_init    = W_PP_pp_end(randperm(N*N));
    W_PPS_po_init   = W_PPS_po_end(randperm(N*N));
    W_PPS_pp_init   = W_PPS_pp_end(randperm(N*N));

    if PP_lesion==1
        disp 'Lesion before training'
        W_PP_po_init    = zeros(N,N);
        W_PP_pp_init    = zeros(N,N);
        lambda_PP       = 0;
    end

    %%% RUN CONSOLIDATION CYCLES
    disp 'running consolidation cycles...'
    % reference memory for consolidation cycles
    % location of object
    p_x    = 0.1768;
    p_y    = -0.1768;
    % object identity
    o_ind  = 1; % this is the object that we will track
    % SC update
    W_SC_po_update  = A_CA3_max*exp(-((p_x-CA3_p_pf(:,1)).^2 + (p_y-CA3_p_pf(:,2)).^2)./(2*CA3_pf_var)) * CA3_o_A(:,o_ind)';
    W_SC_po_update  = W_SC_po_update./max(sum(W_SC_po_update,2));
    W_SC_po_init    = W_SC_po_init + W_SC_po_update*SC_update_fac/(1-SC_update_fac);
    W_SC_po_init    = W_SC_po_init./max(sum(W_SC_po_init,2));

    file = sprintf('%sdata_nsec_%d_seed_%d_cycle_0',dir_res,nsec,seed);
    save(file,'N','CA3_o_A','CA3_p_pf','CA3_pf_var','A_CA3_max','EC_o_A','EC_p_gf','A_EC_max','W_SC_po_init','W_PP_po_init','W_PPS_po_init');

    %% run simulation
    for i_cycle=1:N_cycle
        i_cycle
        params = [nsec dt N W_PP_max W_PPS_max W_SC_id W_SCS_id it_del_CA3 it_del_CA1 lambda_PP lambda_PPS alpha tau_LTP tau_LTD env_len CA3_pf_var A_EC_max A_CA3_max N_obj seed];
        [W_PP_po_end W_PPS_po_end W_PP_pp_end W_PPS_pp_end] = integrate_eqns(W_SC_po_init,W_PP_po_init,W_PP_pp_init,W_PPS_po_init,W_PPS_pp_init,params,CA3_p_pf,CA3_o_A,EC_p_gf,EC_o_A);

        file = sprintf('%sdata_nsec_%d_seed_%d_cycle_%d',dir_res,nsec,seed,i_cycle);
        save(file,'W_SC_po_init','W_PP_po_end','W_PPS_po_end');

        %%% UPDATE INITIAL WEIGHT MATRICES
        % location of object in circular environment
        p_ang  = rand*2*pi;
        p_dis  = rand + rand;
        if p_dis>1, p_dis = 2-p_dis; end
        p_x  = env_len/2*p_dis*cos(p_ang);
        p_y  = env_len/2*p_dis*sin(p_ang);
        % object identity in square environment
        o_ind  = randi([2 N_obj]); % to avoid to hit the memory that we will track (no point here)
        % SC update
        A_CA1_p_tmp     = A_CA3_max*exp(-((p_x-CA3_p_pf(:,1)).^2 + (p_y-CA3_p_pf(:,2)).^2)./(2*CA3_pf_var));
        W_SC_po_update  = A_CA1_p_tmp(randperm(N)) * CA3_o_A(:,o_ind)'; % REMAP CA1_p PATTERN
        W_SC_po_update  = W_SC_po_update./max(sum(W_SC_po_update,2));
        W_SC_po_init    = W_SC_po_init + W_SC_po_update*SC_update_fac/(1-SC_update_fac);
        W_SC_po_init    = W_SC_po_init./max(sum(W_SC_po_init,2));
        W_PP_po_init    = W_PP_po_end;
        W_PP_pp_init    = W_PP_pp_end;
        W_PPS_po_init   = W_PPS_po_end;
        W_PPS_pp_init   = W_PPS_pp_end;

        if PP_lesion==2 && i_cycle>=21
            disp 'Lesion after cycle 21'
            W_PP_po_init    = zeros(N,N);
            W_PP_pp_init    = zeros(N,N);
            lambda_PP       = 0;
        end
    end
end

