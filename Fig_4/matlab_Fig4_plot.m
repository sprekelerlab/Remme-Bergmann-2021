% Memory consolidation - Figure 4
% Michiel Remme, May 2016
% memory consolidation in a hierarchical network
% analyze and plot data after running simulation with matlab_Fig_4_run.m

clear all
clf

%%
seed_init   = 100;
N_cycle     = 100; % (2000) simulation in articles uses 1050 cycles (2000 in total, first 950 are initialization)
cycle_start = 1; %

N_layer     = 8; % (8)
N_cell      = 256;

w_max       = 2/256;

corr_mref_all   = zeros(N_cycle,N_layer+1); % includes the HPC

ind = 1;
for cycle_i = cycle_start+(0:N_cycle-1)
    file    = sprintf('_results/data_Fig4_seed_%d_Nlayer_%d_cycle_%d',seed_init,N_layer,cycle_i);
    if cycle_i==cycle_start
        load(file,'A_','W_HPC_mem','W_HPC_init','W_')
        M_ref       = W_HPC_mem; % this is the reference memory that will be tracked
        W_HPC_ref   = W_HPC_init;
    else
        load(file,'A_','W_','W_HPC_init');
    end
    corr_mref_all(ind,1) = corr(W_HPC_init(:),M_ref(:)); % correlation of HPC matrix with reference memory
    for k = 1:N_layer
        w_tmp = W_(:,:,k);
        corr_mref_all(ind,k+1) = corr(w_tmp(:),M_ref(:)); % correlation of shortcut matrices with reference memory
    end
    ind = ind + 1;
end
t_cycles = 1:size(corr_mref_all,1);

%%

subplot(311)
loglog(t_cycles,corr_mref_all)
axis([1 500 0.01 1])
xlabel('Consolidation cycles')
ylabel('Correlations')
set(gca,'XTick',[1 10 100 1000],'XTickLabel',[1;10;100;1000])
set(gca,'YTick',[0.001 0.01 0.1 1],'YTickLabel',[0.001;0.01;0.1;1])


