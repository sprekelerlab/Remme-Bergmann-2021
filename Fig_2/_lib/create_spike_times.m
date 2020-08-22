%% create place-field sc input and grid-field pp input spike times
% Michiel Remme, May 2016

tend = 5*60; % sec - create spiketimes for tend seconds
Nind = tend*1000/input_vec_dt;

%%%% pp output
Rate_pp = 10; % Hz : mean poisson rate when active

pp_spacing_vec = linspace(2,6,Ne_pp); % how many grid fields per cycle
pp_phase = zeros(Ne_pp,1); % phase shift of grid fields of pp fibers

for k = 1:Ne_pp
    pp_field_size = T_cycle*1000/input_vec_dt/(2*pp_spacing_vec(k));
    pp_phase(k) = randi(floor( 2*pp_field_size )) - pp_field_size; % also have negative phase shifts
    
    gf_vec = zeros(T_cycle*1000/input_vec_dt,1);
    gf_vec_ind = 1:(length(gf_vec)+pp_field_size);
    tmp = mod(ceil((gf_vec_ind-pp_phase(k))/pp_field_size),2);
    gf_vec = tmp(1:length(gf_vec));
    
    pp_sp_tmp = unique(randi(Nind,poissrnd(Rate_pp*tend),1)); % indices van spiketimes als altijd in grid field
    % remove spike times that are outside of grid field
    ind_tmp = mod(pp_sp_tmp,T_cycle*1000/input_vec_dt) + 1; % spike time index within one cycle
    pp_sp_ind{k} = pp_sp_tmp(gf_vec(ind_tmp)==1); % only take those spikes that are in active grid field
end

% now reverse spike order in every 2nd cycle
for k = 1:Ne_pp
    pp_sp_ind_2nd_cycle = mod(pp_sp_ind{k},2*T_cycle*1000/input_vec_dt)>T_cycle*1000/input_vec_dt; % logical vector for all spikes times in 2nd cycles
    pp_sp_ind{k}(pp_sp_ind_2nd_cycle) = (2*ceil(pp_sp_ind{k}(pp_sp_ind_2nd_cycle)/(T_cycle*1000/input_vec_dt)) - 1)*(T_cycle*1000/input_vec_dt) - pp_sp_ind{k}(pp_sp_ind_2nd_cycle);
end

% sort spiketimes, add sentinel for Neuron scanuntil function and turn into spike times instead of indices
sentinel = -1;
for k = 1:Ne_pp
    pp_sp_ind{k} = [sort(pp_sp_ind{k})*input_vec_dt; sentinel];
end

% turn pp indices into arrays that are useful for Neuron
pp_fibs = vertcat(pp_sp_ind{:}); % all spike times concatenated in one vector


%%%% sc output
Rate_sc = 10; % Hz

sc_spacing_vec = 5*ones(Ne_sc,1); % relative size of place field is 1/(2*sc_spacing_vec)
sc_field_size =  T_cycle*1000/input_vec_dt/(2*sc_spacing_vec(1));
sc_phase = zeros(Ne_sc,1); % phase shift of place fields of sc fibers
sc_phase = randi(floor(T_cycle*1000/input_vec_dt + sc_field_size),Ne_sc,1) - sc_field_size; % also have negative phase shifts

for k = 1:Ne_sc
    pf_vec = zeros(T_cycle*1000/input_vec_dt,1);
    pf_vec(max(1,sc_phase(k)):min(sc_phase(k)+sc_field_size,T_cycle*1000/input_vec_dt)) = 1;
    sc_sp_tmp = unique(randi(Nind,poissrnd(Rate_sc*(tend+1)),1)); % indices van spiketimes als altijd in place field (add 1 sec because is subtracted at end
    % remove spike times that are outside of place field
    ind_tmp = mod(sc_sp_tmp,T_cycle*1000/input_vec_dt) + 1; % spike time index within one cycle
    sc_sp_ind{k} = sc_sp_tmp(pf_vec(ind_tmp)==1); % only take those spikes that are in active grid field
end

% now reverse spike order in every 2nd cycle
for k = 1:Ne_sc
    sc_sp_ind_2nd_cycle = mod(sc_sp_ind{k},2*T_cycle*1000/input_vec_dt)>T_cycle*1000/input_vec_dt; % logical vector for all spikes times in 2nd cycles
    sc_sp_ind{k}(sc_sp_ind_2nd_cycle) = (2*ceil(sc_sp_ind{k}(sc_sp_ind_2nd_cycle)/(T_cycle*1000/input_vec_dt)) - 1)*(T_cycle*1000/input_vec_dt) - sc_sp_ind{k}(sc_sp_ind_2nd_cycle);
end

% sort spiketimes, add sentinel for NEURON scanuntil function and turn into spike times instead of indices
sentinel = -1;
for k = 1:Ne_sc
    sc_sp_ind{k} = [sort(sc_sp_ind{k})*input_vec_dt; sentinel];
end

% turn sc indices into arrays that are useful for Neuron
sc_fibs = vertcat(sc_sp_ind{:}); % all spike times concatenated in one vector


%%%% save initialization files
save('../_initialization_files/place_grid_fields_2500_500', 'pp_phase', 'pp_spacing_vec', 'T_cycle', 'sc_phase', 'sc_spacing_vec','input_vec_dt');

filename = '../_initialization_files/spiketimes_2500_500_5min_pp.dat';
fid = fopen(filename,'w');
fprintf(fid,'%.2f\n',pp_fibs);
fclose(fid);

filename = '../_initialization_files/spiketimes_2500_500_5min_sc.dat';
fid = fopen(filename,'w');
fprintf(fid,'%.2f\n',sc_fibs);
fclose(fid);

