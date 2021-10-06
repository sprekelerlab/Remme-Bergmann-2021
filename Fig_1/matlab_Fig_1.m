% Memory consolidation - Figure 1
% Michiel Remme, May 2016
% Single compartment model

clear all

% the gnu scientific library (gsl) needs to be installed because of the random number generator used in the c-code

% compile the c-code
% on mac
mex -I/usr/local/include -L/usr/local/lib -lgsl -lm integrate_eqns.c
% on linux machines need to add lgslcblas
% mex -lgsl -lgslcblas -lm integrate_eqns.c

% octave: mkoctfile --mex -I/usr/local/include -L/usr/local/lib -lgsl -lm integrate_eqns.c

%% Parameters

nsec        = 300*60;                   % sec - total duration of simulation (in Fig1: 300*60 sec)
dt          = 0.1;                      % msec

Ne          = 1000;                     % number of excitatory inputs per pathway
E_freq      = 10;                       % Hz - mean input rate

pp_plastic  = true;                     % if not true SC is plastic

% learning parameters
delay_sc_pp = 5;                        % ms - delay between SC and PP pathway
we_max      = 0.006;                    % dimensionless max weight of inputs (i.e., it's relative to membrane leak conductance)
lambda_pp   = pp_plastic*0.005*we_max;  % learning rate PP inputs
lambda_sc   = (1-pp_plastic)*0.005*we_max; % learning rate SC inputs
alpha       = 1.05;                     % ratio of LTD vs LTP

save_we_ivl = 30;                       % sec - save weight matrices every X seconds

we_mu       = 0.05;                     % mu-parameter for the mixture of exponentials weight distribution

% set intial weights of inputs
if pp_plastic
    we_tmp      = bimodal_dist(we_mu, Ne/2); % create vector of random weights with bimodal distribution
    we_sc_init  = [we_max*we_tmp(randperm(Ne/2)); we_max*we_tmp(Ne/2+randperm(Ne/2))];
    we_tmp      = bimodal_dist(we_mu, Ne/2);
    we_pp_init  = we_max*we_tmp(randperm(Ne));
else
    we_tmp      = bimodal_dist(we_mu, Ne/2);
    we_sc_init  = we_max*we_tmp(randperm(Ne));
    we_tmp      = bimodal_dist(we_mu, Ne/2);
    we_pp_init  = [we_max*we_tmp(randperm(Ne/2)); we_max*we_tmp(Ne/2+randperm(Ne/2))];
end

%% Run c-code
params = [nsec dt Ne E_freq we_max delay_sc_pp lambda_sc lambda_pp alpha save_we_ivl];
[we_sc_ we_pp_ rate_ v_] = integrate_eqns(we_sc_init,we_pp_init,params);

%% Plot results

t  = (1:nsec)/60;
t_ = (save_we_ivl:save_we_ivl:nsec)/60;

we_sc_end = we_sc_(:,end);
we_pp_end = we_pp_(:,end);

% compute correlations bewteen SC and PP pathway
corr_pp_sc = [];

for ind = 1:nsec/save_we_ivl
    wsc_tmp      = we_sc_(:,ind) - mean(we_sc_(:,ind));
    wpp_tmp      = we_pp_(:,ind) - mean(we_sc_(:,ind));
    corr_pp_sc(ind) = (wpp_tmp'*wsc_tmp) / (sqrt(wpp_tmp'*wpp_tmp)*sqrt(wsc_tmp'*wsc_tmp));
end


figure(1)

subplot(122)
plot(t_,corr_pp_sc,'k')
axis([0 nsec/60 -0.2 1.01])
xlabel('Time (min)')
ylabel('Correlation')

subplot(221)
plot(1:Ne,we_pp_end./we_max,'r.')
xlabel('PP inputs')
ylabel('Normalized weights')
set(gca,'xtick',[0:200:Ne])
set(gca,'ytick',[0:0.2:1])
axis([0 Ne 0 1])

subplot(223)
plot(1:Ne,we_sc_end./we_max,'b.')
xlabel('SC inputs')
ylabel('Normalized weights')
set(gca,'xtick',[0:200:Ne])
set(gca,'ytick',[0:0.2:1])
axis([0 Ne 0 1])

% pause()
