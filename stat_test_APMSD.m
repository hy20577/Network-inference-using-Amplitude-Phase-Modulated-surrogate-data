function [TPR_mat_heatmap, FPR_mat_heatmap, inf_nets, Mir_overall, output_all] = stat_test_APMSD(Orbits, Adj, p_sl, varargin)

% This function implements the APMSD method for the given time recording. 

% INPUTs: 

% Orbits: Matrix of all dataset, N-by-M for N time recording and M nodes.

% Adj : Adjacency matrix of original network

% p_sl: significance level of hypothesis tests, default 0.1.

% 'varargin' can take 2 values for pc_{1} and pc_{2}. Its optional values
% 0:0.1:1. This produces 11-by-11 parameter space in the paper.

% OUTPUTs


if ~exist("p_sl", "var")
p_sl = 0.1;
end

if nargin == 2
    pc1 = 0:0.1:1;
    pc2 = 0:0.1:1;
elseif nargin == 3
    pc1 = 0:0.1:1;
    pc2 = 0:0.1:1;
elseif nargin == 4
    pc1 =  varargin{1};
    pc2 = 0:0.1:1;
elseif nargin == 5
    pc1 = varargin{1};
    pc2 = varargin{2};
end

time_range = size(Orbits,1);
M = size(Adj,2);

l1 = length(pc1);
l2 = length(pc2);
output_all = zeros(l1*l2, 5); % pc1, pc2, TPR, FPR, distance  
TPR_matrix = zeros(l1,l2);
FPR_matrix = zeros(l1,l2);
dens = zeros(l1,l2);
dens_orig = sum(sum(Adj))/(M*(M-1));
counter =0;

inf_nets = zeros(M,M,l1*l2);
Mir_overall = zeros(1/p_sl+1,M*(M-1)/2,l1*l2);

for i=1:l1
    for j=1:l2
        counter = counter+1;
        [MIR_all, adj_MIR_fdr, TP_FP_Rate] = Mir_surrogate_data(Orbits, p_sl, time_range, "APMSD", Adj,'', pc1(i), pc2(j));
        Mir_overall(:,:,counter) = MIR_all;
        output_all(counter,1) = pc1(i);
        output_all(counter,2) = pc2(j);
        output_all(counter,3:4) = TP_FP_Rate(2:3);
        output_all(counter,5) = sqrt((1-TP_FP_Rate(2))^2+TP_FP_Rate(3)^2);

        TPR_matrix(i,j) = TP_FP_Rate(2);
        FPR_matrix(i,j) = TP_FP_Rate(3);
        dens(i,j) = sum(sum(adj_MIR_fdr))/(M*(M-1)) - dens_orig; % density of the inferred network - original network.
        
        inf_nets(:,:,counter) = adj_MIR_fdr;

        sprintf('%d out of %d has done', counter, l1*l2)
    end
end
       output_all = sortrows(output_all,5);

TPR_mat_heatmap = zeros(l1,l2);
FPR_mat_heatmap = zeros(l1,l2);

for i=1:l1
TPR_mat_heatmap(i,:) = TPR_matrix(l1-i+1, :); % upside down TPR and FPR to present results on y-axis in descending order.
FPR_mat_heatmap(i,:) = FPR_matrix(l1-i+1, :);
end

end