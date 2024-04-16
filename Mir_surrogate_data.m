
function [MIR_overall, adj_MIR_fdr, TP_FP_Rate] = Mir_surrogate_data(All_Dat, p_sl, time_range, type_of_surr, original_network, label, varargin)

% This function generates surrogate data, compute their MIR values and
% compare them with the original MIR values to obtain probability
% statistics for all pairs in the network. After false discovery rate (FDR)
% process, it produces the final inferred network adj_MIR_fdr. 


%%% COMPULSORY INPUTs

% All_Dat: Data to use computation of MIR after removing transient period.

% type_of_surr: "twnsd", "randsd", "APMSD"

%%% OPTIONAL INPUTs 

% p_sl: confidence interval for hypothesis testing (default: 0.01).

% desired_corr: When the surrogate data is generated as correlated between
        % pairs, you need to define desired correlation otherwise it is
        % defined as 0.5 as default.

% time_range : To consider different time length. Default is
% 'logspace(3e2,size(All_Dat,1),100)';


% OUTPUTs

    % MIR_overall : First row gives MIR value for original data. Next
           % 1/p_sl row gives the MIR values for surrogate data.

    % adj_MIR_fdr : gives the adjacency matrix (symetric) after fdr
                % process.

    % TP_FP_Rate : nx3 matrix. 1st column -- Time, 2nd -- TPR, 3rd -- FPR.
    
    
time_start = tic;

if ~exist('label', 'var') % If label does not exist, assign empty string.
label = '';
end

% if nargin == 6  
%      pc1 = 1;
%      pc2 = 0.5;
% elseif nargin == 7
%    pc1 = varargin{1};
%    pc2 = 0.5;
% elseif nargin ==8
%     pc1 = varargin{1};
%     pc2 = varargin{2};
% end

% To consider different time length in twin surrogate and random surrogate
% methods. We consider whole dataset for APMSD method. 

if ~exist("time_range",'var')
time_range = logspace(3e2,size(All_Dat,1),100);
end

if time_range(end) > size(All_Dat,1)
    time_range(end) = size(All_Dat,1);
    warning('the last element of time_range should be less than number of observations in original data, assigned the biggest possible value automatically')
end

if ~exist("p_sl",'var') % Significance level of hypothesis tests, default value 0.1.
    p_sl = 0.1;
end

if type_of_surr == "APMSD"  % if pc_{1} and pc_{2} parameters are not defined in APMSD method, but we always define it stat_test_APMSD function
         if nargin == 6
             pc1 = 1;
             pc2 = 0.5;
         elseif nargin == 7
           pc1 = varargin{1};
           pc2 = 0.5;
         elseif nargin ==8
            pc1 = varargin{1};
            pc2 = varargin{2};
         end
end

nVar = size(All_Dat,2);
Nsd = 1/p_sl;   % Number of surrogate data
MIR_overall = zeros(1+Nsd,nVar*(nVar-1)/2, length(time_range)); % initiate the output matrix.
Inferred_Adj_MIR = zeros(nVar, nVar,length(time_range)); % initiate the matrix.
adj_MIR_fdr = zeros(nVar,nVar, length(time_range)); % initiate the output matrix.

TP_FP_Rate = zeros(length(time_range),3);
TP_FP_Rate(:,1) = time_range;

Dur = zeros(length(time_range),1); % initiate the time length vector


for ITER=1:length(time_range)

    Dat = All_Dat(1:time_range(ITER),:);

    [MIR_orig, ~] = MIR(Dat);

    niter = 0;

    for n = 1:nVar
        for m = n+1:nVar
          niter = niter + 1;
          MIR_overall(1,niter,ITER) = MIR_orig(n,m);
        end
    end

    % Producing normal surrogate data and their MI values.
MIR_surrogate = zeros(Nsd, nVar*(nVar-1)/2);

parfor i=1:Nsd  % Generate surrogate data and calculate their MIR values.
     
     surr_dat = zeros(size(Dat));   
     if type_of_surr == "twnsd"
        surr_dat = 	phaseran(Dat,1);  % phaserun produces twin surrogate data to test the null hypotheses      
     elseif  type_of_surr == "randsd"
        surr_dat = rand(length(Dat), nVar); % generate random surrogate data
     elseif type_of_surr == "APMSD"         
            surr_dat = APMSD(Dat, 1, pc1, pc2); % generate APMSD
     end

       [MIR_surr, ~] = MIR(surr_dat);
       niter=0;
       MIR_array = zeros(1,nVar*(nVar-1)/2);

       for n = 1:nVar
            for m = n+1:nVar
              niter = niter+1;
              MIR_array(niter) = MIR_surr(n,m);
            end
       end  
        MIR_surrogate(i,:) = MIR_array;
end

    MIR_overall(2:end,:,ITER) = MIR_surrogate;
    c_xy_MIR = zeros(1, nVar*(nVar-1)/2); % number of surrogate data whose MIRs are equal to or higher than original MIR.
    
    for i=1:nVar*(nVar-1)/2
        sum1 = 0;
       for j=2:Nsd+1

          if MIR_overall(j,i,ITER) > MIR_overall(1,i,ITER) || MIR_overall(j,i,ITER) == MIR_overall(1,i,ITER)
              sum1 = sum1 + 1;
          end
          c_xy_MIR(i) = sum1;
       end
    end

    prob_MIR = c_xy_MIR / Nsd;

    iter = 0;

    for i=1:nVar
         for j= i+1:nVar
          iter = iter+1;

          if prob_MIR(iter) < p_sl
                Inferred_Adj_MIR(i,j,ITER) = 1;
          end

         end
    end

    Inferred_Adj_MIR(:,:,ITER) = Inferred_Adj_MIR(:,:,ITER) + Inferred_Adj_MIR(:,:,ITER)';

    % Calculation of False Discovery Rate (FDR)

     [adj, ~, ~, ~]=fdr_bh(prob_MIR,p_sl/100,'pdep','no');
     
     counter = 0;
     for i=1:nVar
        for j=i+1:nVar
          counter = counter+1;
          adj_MIR_fdr(i,j,ITER) = adj(counter);
        end
     end
  
     adj_MIR_fdr(:,:,ITER) = adj_MIR_fdr(:,:,ITER) + adj_MIR_fdr(:,:,ITER)'; 
    [TP_Rate,FP_Rate] = TPR_FPR(original_network,adj_MIR_fdr(:,:,ITER));
    TP_FP_Rate(ITER,2) = TP_Rate;
    TP_FP_Rate(ITER,3) = FP_Rate; 

    Dur(ITER) = toc(time_start);


%%% Let's consider the data length in each iteration.

rate_of_completed = sum(time_range(1:ITER))/sum(time_range);
remaining_time_seconds = ((Dur(ITER))/rate_of_completed)-(Dur(ITER));
 
fprintf("Completed: %s %% \n ETA: %s (in DD:HH:MM:SS.MS)\n", num2str(round(rate_of_completed*100,1)), datestr(remaining_time_seconds/(24*60*60),'DD:HH:MM:SS.FFF') )


end

fprintf('Elapsed time (function:Mir_surrogate_data) for %s and p_sl: %g = %s (in DD:HH:MM:SS.MS)\n',label,p_sl, datestr(toc(time_start)/(24*60*60),'DD:HH:MM:SS.FFF'))

end
