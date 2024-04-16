
function [Mir, Mi] = MIR(Dat)

%%% The function computes MIR values (without normalisation as in Bianco et.
% al. 2016) and chose the maximum partition size N which satisfies the 
% unbiased probability calculation in eq.17 in the same paper. 

%%% INPUT:

% Dat: N-by-M matrix for N time recording and M variables. 

%%% OUTPUT

% Mir : Mutual Information Rate

% Mi : Mutual Information

M = length(Dat(1,:));
[N_min, N_max] = Grid(Dat);

Is = zeros(M,M);  %% Mutual Information between pairs
Ic = zeros(M,M);  %% MIR between pairs gaining by dividing MI by t.
Mi_MIR_over_partitions = zeros(2,length(N_min:N_max));

for N = N_min:N_max
    
I_s = zeros(M,M);  %% Mutual Information between pairs
I_c = zeros(M,M);  %% MIR between pairs gaining by dividing MI by t.
    
      for i=1:M-1  %par
        for j=(i+1):M
        node1 = Dat(:,i);
        node2 = Dat(:,j);
        node12 = [node1, node2];
        [coordinates, ppb] = location(node12,N);        
        t = corr_decay(node12, N); %% correlation decay time for each pair.
        t = t(t~=0);
       
        MI=0;
               for m=1:N^2
               prob_node1 = sum(coordinates(:,2) == ceil(m/N)) / length(node1);
                 if mod(m,N) ~= 0
                    prob_node2 = sum(coordinates(:,3) == mod(m,N))/length(node1);
                 else 
                    prob_node2 = sum(coordinates(:,3) == N)/length(node1);
                 end
               joint_prob = ppb(m,2)/length(node1);
                         if joint_prob ~= 0  %%  p*log(p) -> 0 assumed as p ->0.
                         MI = MI + joint_prob * log(joint_prob/(prob_node1*prob_node2));
                         end             
               end
       
             I_s(i,j) =  MI;
             I_c(i,j) =  MI/t;
             
        end
      end
 Is = Is + I_s ;  
 Ic = Ic + I_c ; 
           
end
      
Mir = Ic/length(N_min:N_max);
Mi = Is/length(N_min:N_max);
 if exist('nodes','var')
     tiledlayout(1,2)
     nexttile
     loglog(N_min:N_max, Mi_MIR_over_partitions(1,:),'o--')
     title('Evolution of MI over partitions')
     nexttile
     loglog(N_min:N_max, Mi_MIR_over_partitions(2,:),'.-')
     title('Evolution of MIR over partitions')
 end

end