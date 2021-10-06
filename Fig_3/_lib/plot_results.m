% Michiel Remme, May 2016
% code to analyze results when inputs are initialized with SC tuned to place field

% load results from NEURON
we_sc_      = load(fullfile('..',dir_res,'wesc_mat.dat'));
we_pp_      = load(fullfile('..',dir_res,'wepp_mat.dat'));
spikes_sc   = load(fullfile('..',dir_res,'spikes_vec_sc.dat'));
spikes_pp   = load(fullfile('..',dir_res,'spikes_vec_pp.dat'));
vdata       = load(fullfile('..',dir_res,'vdata.dat'));

vs      = vdata(:,1); % soma V
vd      = vdata(:,2); % dend V
t       = vdata(:,3);

nsec    = tstop/1000;
t_sec   = 0:nsec;
t_min   = (1:nsec)/60; % in minutes
t_mat   = (0:save_we_ivl:nsec)/60; % times that weights have been saved

rate_sc_ = histc(spikes_sc/1000,t_sec);
rate_pp_ = histc(spikes_pp/1000,t_sec);

we_sc_init  = we_sc_(:,1);
we_pp_init  = we_pp_(:,1);
we_sc_end   = we_sc_(:,end);
we_pp_end   = we_pp_(:,end);

%% SC INPUTS
% weighted sum of tuning functions
sc_tuning = zeros(1000,size(we_sc_,2)+1);
sc_tuning_tmp = zeros(1000,1);
for k = 1:Ne_sc
    ind_vec_start = max(1,round(sc_phase(k)./(T_cycle/input_vec_dt) + 1));
    ind_vec_end = min(round(sc_phase(k)./(T_cycle/input_vec_dt) + 1000./(sc_spacing_vec(k)*2) + 1),length(sc_tuning_tmp));
    sc_tuning_tmp(ind_vec_start:ind_vec_end) = sc_tuning_tmp(ind_vec_start:ind_vec_end) + we_sc_init(k);
end
sc_tuning(:,1) = sc_tuning_tmp(1:1000);
for kk = 2:size(sc_tuning,2)
    we_kk = kk-1;
    sc_tuning_tmp = zeros(2000,1);
    for k = 1:Ne_sc
        ind_vec_start = max(1,round(sc_phase(k)./(T_cycle/input_vec_dt) + 1));
        ind_vec_end = min(round(sc_phase(k)./(T_cycle/input_vec_dt) + 1000./(sc_spacing_vec(k)*2) + 1),length(sc_tuning_tmp));
        sc_tuning_tmp(ind_vec_start:ind_vec_end) = sc_tuning_tmp(ind_vec_start:ind_vec_end) + we_sc_(k,we_kk);
    end
    sc_tuning(:,kk) = sc_tuning_tmp(1:1000);
end

%% PP INPUTS
% weighted sum of tuning functions
pp_tuning = zeros(1000,size(we_pp_,2)+1);
pp_tuning_tmp = zeros(1000,1);
for k = 1:Ne_pp
    pp_field_size = 1000/(2*pp_spacing_vec(k));
    gf_vec = zeros(1000,1);
    gf_vec_ind = 1:(length(gf_vec)+pp_field_size);
    tmp = mod(ceil((gf_vec_ind-pp_phase(k)/(T_cycle/input_vec_dt))/pp_field_size),2);
    gf_vec = we_pp_init(k)*tmp(1:length(gf_vec))';
    pp_tuning_tmp = pp_tuning_tmp + gf_vec;
end
pp_tuning(:,1) = pp_tuning_tmp;
for kk = 2:size(pp_tuning,2)
    we_kk = kk-1;
    pp_tuning_tmp = zeros(1000,1);
    for k = 1:Ne_pp
        pp_field_size = 1000/(2*pp_spacing_vec(k));
        gf_vec = zeros(1000,1);
        gf_vec_ind = 1:(length(gf_vec)+pp_field_size);
        tmp = mod(ceil((gf_vec_ind-pp_phase(k)/(T_cycle/input_vec_dt))/pp_field_size),2);
        gf_vec = we_pp_(k,we_kk)*tmp(1:length(gf_vec))';
        pp_tuning_tmp = pp_tuning_tmp + gf_vec;
    end
    pp_tuning(:,kk) = pp_tuning_tmp;
end

%% INPUT TUNING

% phase of ca1 'place field'
ca1_phase   = 0.5*T_cycle*1000/input_vec_dt; % time in cycle
pp_grid_size = (T_cycle*1000/input_vec_dt)./(pp_spacing_vec*2);
pp_pot_vec  = mod(floor((ca1_phase-pp_phase')./pp_grid_size),2)==0;

% PP FIBERS
% weighted sum of tuning functions
pp_tuning_theor = zeros(1000,1);
pp_tuning_tmp = zeros(1000,1);
for k = 1:Ne_pp
    % 'optimal' grid field
    pp_field_size = 1000/(2*pp_spacing_vec(k));
    gf_vec = zeros(1000,1);
    gf_vec_ind = 1:(length(gf_vec)+pp_field_size);
    tmp = mod(ceil((gf_vec_ind-pp_phase(k)/(T_cycle/input_vec_dt))/pp_field_size),2);
    gf_vec = tmp(1:length(gf_vec))';

    pp_tuning_tmp = pp_tuning_tmp + pp_pot_vec(k)*we_pp_max*gf_vec;
end
pp_tuning_theor(:,1) = pp_tuning_tmp;
pp_pot_end = we_pp_end'>we_pp_max/2; % vector of potentiated inputs at end of simulation

% compute correlation between we_pp and 'ideal pp weight vector'
we_pp_corr = zeros(size(we_pp_,2),1);
we_pp_ideal = pp_pot_vec'*we_pp_max;
for k = 1:size(we_pp_,2)
    we_tuning_pearson(k) = corr(pp_tuning(:,k),pp_tuning_theor); % pearson correlation coeff (with means removed)
end

%%

figure(1)
clf

subplot(311)
plot((0:999)*1e-3,sc_tuning(:,1)*1e3,'b',(0:999)*1e-3,pp_tuning(:,1)*1e3,'r')
axis([0 1 0 90])
xlabel('Position')
ylabel('Input tuning (nS)')
set(gca,'xtick',[0:0.2:1])
set(gca,'ytick',[0:20:120])
title('Before consolidation')

subplot(312)
plot((0:999)*1e-3,sc_tuning(:,end)*1e3,'b',(0:999)*1e-3,pp_tuning_theor*1e3,'m',(0:999)*1e-3,pp_tuning(:,end)*1e3,'r')
axis([0 1 0 90])
xlabel('Position')
ylabel('Input tuning (nS)')
set(gca,'xtick',[0:0.2:1])
set(gca,'ytick',[0:20:120])
title('After consolidation')

subplot(313)
plot(t_mat,we_tuning_pearson,'k')
axis([0 5 -0.05 1])
xlabel('Time (min)')
ylabel('Correlation')
set(gca,'xtick',[0:1:5])
set(gca,'ytick',[0:0.5:1])


