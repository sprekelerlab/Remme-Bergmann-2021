% Memory consolidation - Figure 3
% Michiel Remme, May 2016
% memory consolidation in a hippocampal network
% analyze and plot data after running simulation with matlab_Fig_3_run.m

clear all
clf

% PARAMETERS
Ntrial      = 1; % 10
Ncycle      = 31; % 31
seed_range  = 0+(1:Ntrial); % 50
Ncycle_ivl  = 7; % 7 (only for plotting purposes)
cycle_range = 1:Ncycle_ivl:Ncycle; % which cycles to show
t_cycle     = 0:Ncycle;

nsec        = 150;

env_len     = 1;
env_rad     = env_len/2;
dx_test     = env_len/21; % spatial resolution to infer probabilities

P_SC_var    = 0.3^2; % variance of decoder
P_PP_var    = 0.3^2;
P_PPS_var   = 0.3^2;

pos_vec_test = -env_len/2:dx_test:env_len/2;
[xvec,yvec] = meshgrid(pos_vec_test);
p_x_test    = xvec(:);
p_y_test    = yvec(:);
p_dist_test = sqrt(p_x_test.^2+p_y_test.^2);
quad_1      = find(p_dist_test<env_rad & p_x_test>0 & p_y_test>0);
quad_2      = find(p_dist_test<env_rad & p_x_test<0 & p_y_test>0);
quad_3      = find(p_dist_test<env_rad & p_x_test<0 & p_y_test<0);
quad_4      = find(p_dist_test<env_rad & p_x_test>0 & p_y_test<0);
n_grid_points = length(p_x_test);

% object location
p_x         = 0.1768;
p_y         = -0.1768;
% object identity
o_ind       = 1;

%%
P_SC_quad_cycle  = zeros(4,Ncycle+1,Ntrial);
P_PP_quad_cycle  = zeros(4,Ncycle+1,Ntrial);
P_PPS_quad_cycle = zeros(4,Ncycle+1,Ntrial);
P_SC_prob_cycle  = zeros(length(pos_vec_test),length(pos_vec_test),Ncycle+1,Ntrial);
P_PP_prob_cycle  = zeros(length(pos_vec_test),length(pos_vec_test),Ncycle+1,Ntrial);
P_PPS_prob_cycle = zeros(length(pos_vec_test),length(pos_vec_test),Ncycle+1,Ntrial);

W_SC_mean = zeros(Ncycle,Ntrial);
W_PP_mean = zeros(Ncycle,Ntrial);
W_PPS_mean = zeros(Ncycle,Ntrial);

figure(1)
clf
for trial_i = 1:Ntrial
    fprintf('Trial No. %d\n',trial_i);
    file_i = seed_range(trial_i);
    for cycle_i = 0:Ncycle
        file    = sprintf('_results/data_nsec_%d_seed_%d_cycle_%d',nsec,file_i,cycle_i);
        if cycle_i == 0
            % store CA1 and SUB activity patterns for different input locations at onset (= same for all cycles)
            load(file,'N','CA3_o_A','CA3_p_pf','CA3_pf_var','A_CA3_max','EC_o_A','EC_p_gf','A_EC_max','W_SC_po_init','W_PP_po_init','W_PPS_po_init');
            % CA1_p activity through CA3_p identity input for locations on grid
            A_CA1_p_CA3_p   = zeros(N,n_grid_points);
            for i = 1:n_grid_points
                A_CA1_p_CA3_p(:,i) = A_CA3_max*exp(-((p_x_test(i)-CA3_p_pf(:,1)).^2 + (p_y_test(i)-CA3_p_pf(:,2)).^2)./(2*CA3_pf_var));
                A_CA1_p_CA3_p(:,i) = A_CA1_p_CA3_p(:,i)-mean(A_CA1_p_CA3_p(:,i)); % subtract mean activity
                A_CA1_p_CA3_p(:,i) = A_CA1_p_CA3_p(:,i)./std(A_CA1_p_CA3_p(:,i)); % normalize by std activity
            end
            % SUB_p activity through CA3_p -> CA1_p identity input for locations on grid
            A_SUB_p_CA3_p   = A_CA1_p_CA3_p;
            A_CA3_o         = CA3_o_A(:,o_ind);         % CA3_o activity for stored object
            A_EC_o          = EC_o_A(:,o_ind);          % EC_o activity for stored object
            num_neurons     = size(A_EC_o, 1);
            W_PP_po         = reshape(W_PP_po_init, [num_neurons, num_neurons]);
            W_PPS_po        = reshape(W_PPS_po_init, [num_neurons, num_neurons]);
        else
            load(file,'W_SC_po_init','W_PP_po_end','W_PPS_po_end');
            W_PP_po     = W_PP_po_end;
            W_PPS_po    = W_PPS_po_end;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ANALYZE W_SC MATRIX

        A_CA1_p_CA3_o   = W_SC_po_init*A_CA3_o;              % CA1_p activity through CA3_o input
        A_CA1_p_CA3_o   = A_CA1_p_CA3_o-mean(A_CA1_p_CA3_o); % subtract mean activity
        A_CA1_p_CA3_o   = A_CA1_p_CA3_o./std(A_CA1_p_CA3_o); % normalize by std activity

        P_SC_po         = zeros(n_grid_points,1); % probability of object at position p
        for j = 1:n_grid_points
            P_SC_po(j)  = exp(-( sum((A_CA1_p_CA3_p(:,j) - A_CA1_p_CA3_o).^2)/N )./(2*P_SC_var));
        end
        P_SC_po         = P_SC_po./sum(P_SC_po); % normalize sum to 1
        P_SC_quad(1)    = sum(P_SC_po(quad_1));
        P_SC_quad(2)    = sum(P_SC_po(quad_2));
        P_SC_quad(3)    = sum(P_SC_po(quad_3));
        P_SC_quad(4)    = sum(P_SC_po(quad_4));
        P_SC_quad       = P_SC_quad./sum(P_SC_quad); % because some probs are outside circle
        P_SC_quad_cycle(:,cycle_i+1,trial_i) = P_SC_quad;
        P_SC_prob_cycle(:,:,cycle_i+1,trial_i) = reshape(P_SC_po,length(pos_vec_test),length(pos_vec_test));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ANALYZE W_PP MATRIX

        A_CA1_p_EC_o    = W_PP_po*A_EC_o;          % CA1_p activity through EC_o input for different moments during consolidation
        A_CA1_p_EC_o    = A_CA1_p_EC_o-mean(A_CA1_p_EC_o);  % subtract mean activity
        A_CA1_p_EC_o    = A_CA1_p_EC_o./std(A_CA1_p_EC_o);  % normalize by std activity

        P_PP_po         = zeros(1,n_grid_points); % probability of object at position p
        P_PP_quad       = zeros(1,4); % probs in quadrants
        for j = 1:n_grid_points
            P_PP_po(j)  = exp(-( sum( (A_CA1_p_CA3_p(:,j) - A_CA1_p_EC_o).^2 )/N )./(2*P_PP_var));
        end
        P_PP_po         = P_PP_po./sum(P_PP_po); % normalize sum to 1
        P_PP_quad(1)    = sum(P_PP_po(quad_1));
        P_PP_quad(2)    = sum(P_PP_po(quad_2));
        P_PP_quad(3)    = sum(P_PP_po(quad_3));
        P_PP_quad(4)    = sum(P_PP_po(quad_4));
        P_PP_quad       = P_PP_quad./repmat(sum(P_PP_quad,2),1,4); % because some probs are outside circle
        P_PP_quad_cycle(:,cycle_i+1,trial_i) = P_PP_quad(end,:);
        P_PP_prob_cycle(:,:,cycle_i+1,trial_i) = reshape(P_PP_po,length(pos_vec_test),length(pos_vec_test));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ANALYZE W_PPS MATRIX

        A_SUB_p_EC_o    = W_PPS_po*A_EC_o;        % SUB_p activity through EC_o input for different moments during consolidation
        A_SUB_p_EC_o    = A_SUB_p_EC_o-mean(A_SUB_p_EC_o); % subtract mean activity
        A_SUB_p_EC_o    = A_SUB_p_EC_o./std(A_SUB_p_EC_o); % normalize by std activity

        P_PPS_po    = zeros(1,n_grid_points); % probability of object at position p
        P_PPS_quad  = zeros(1,4); % probs in quadrants
        for j = 1:n_grid_points
            P_PPS_po(j) = exp(-( sum( (A_SUB_p_CA3_p(:,j) - A_SUB_p_EC_o).^2 )/N )./(2*P_PPS_var));
        end
        P_PPS_po        = P_PPS_po/sum(P_PPS_po); % normalize sum to 1
        P_PPS_quad(1)   = sum(P_PPS_po(quad_1));
        P_PPS_quad(2)   = sum(P_PPS_po(quad_2));
        P_PPS_quad(3)   = sum(P_PPS_po(quad_3));
        P_PPS_quad(4)   = sum(P_PPS_po(quad_4));
        P_PPS_quad      = P_PPS_quad./repmat(sum(P_PPS_quad,2),1,4); % because some probs are outside circle
        P_PPS_quad_cycle(:,cycle_i+1,trial_i) = P_PPS_quad(end,:);
        P_PPS_prob_cycle(:,:,cycle_i+1,trial_i) = reshape(P_PPS_po,length(pos_vec_test),length(pos_vec_test));
    end

    figure(1)
    subplot(311)
    plot(t_cycle,P_SC_quad_cycle([1 2 3],:,trial_i),'k',t_cycle,P_SC_quad_cycle(4,:,trial_i),'r');
    ylabel('SC')
    axis([0 Ncycle-1 0 1])
    set(gca,'xtick',[0:Ncycle_ivl:Ncycle],'ytick',[0:0.2:1])
    title('Inferred probability of object in quadrant (single trial)')
    box off

    subplot(312)
    plot(t_cycle,P_PP_quad_cycle([1 2 3],:,trial_i),'k',t_cycle,P_PP_quad_cycle(4,:,trial_i),'r');
    ylabel('PP')
    axis([0 Ncycle-1 0 1])
    set(gca,'xtick',[0:Ncycle_ivl:Ncycle],'ytick',[0:0.2:1])
    box off

    subplot(313)
    plot(t_cycle,P_PPS_quad_cycle([1 2 3],:,trial_i),'k',t_cycle,P_PPS_quad_cycle(4,:,trial_i),'r');
    ylabel('PP SUB')
    xlabel('Days after storage')
    axis([0 Ncycle-1 0 1])
    set(gca,'xtick',[0:Ncycle_ivl:Ncycle],'ytick',[0:0.2:1])
    box off
    
    pause(1)
end

%%

figure(2)
clf

example_trial_idx = 1; % which trial to use as an example

% scale all probability maps with same factor for plotting purposes
zq_max              = max([max(P_SC_prob_cycle(:)) max(P_PP_prob_cycle(:)) max(P_PPS_prob_cycle(:))]);
P_SC_prob_cycle     = P_SC_prob_cycle./zq_max;
P_PP_prob_cycle     = P_PP_prob_cycle./zq_max;
P_PPS_prob_cycle    = P_PPS_prob_cycle./zq_max;

color_max = 1;

for i = 1:length(cycle_range)
    subplot(3,length(cycle_range),i)
    imagesc(pos_vec_test,pos_vec_test,P_SC_prob_cycle(:,:,cycle_range(i),example_trial_idx),[0 color_max])
    axis xy
    axis square
    axis([-0.5 0.5 -0.5 0.5])
    axis off

    subplot(3,length(cycle_range),length(cycle_range)+i)
    imagesc(pos_vec_test,pos_vec_test,P_PP_prob_cycle(:,:,cycle_range(i),example_trial_idx),[0 color_max])
    axis xy
    axis square
    axis([-0.5 0.5 -0.5 0.5])
    axis off

    subplot(3,length(cycle_range),2*length(cycle_range)+i)
    imagesc(pos_vec_test,pos_vec_test,P_PPS_prob_cycle(:,:,cycle_range(i),example_trial_idx),[0 color_max])
    axis xy
    axis square
    axis([-0.5 0.5 -0.5 0.5])
    axis off
end
