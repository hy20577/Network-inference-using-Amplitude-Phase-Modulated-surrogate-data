
% This script is designed to show to how to use APMSD method on the data generated from the given adjacency matrix.

% Load data and adjacency matrix.

Dat = readmatrix('Orbits_circle.txt');
M = size(Dat,2); % number of columns
Adj = readmatrix("adj_16.txt");


[TPR_mat_heatmap, FPR_mat_heatmap, inf_nets, Mir_overall, output_all] = stat_test_APMSD(Dat, Adj);

%%

compactheatmap_fromTFPR(TPR_mat_heatmap, FPR_mat_heatmap)