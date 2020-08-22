% Memory consolidation - Figure 2
% Michiel Remme, May 2016
% for plotting results from a previously run simulation

clear all

dir_init = '_initialization_files';
dir_lib  = '_lib';
dir_res  = '_results/';

run(fullfile(dir_init,'params.m'));
load(fullfile(dir_init,'place_grid_fields_2500_500'))% load place/grid field parameters
run(fullfile(dir_lib,'plot_results'))
