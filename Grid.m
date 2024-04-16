function [N_min, N_max] = Grid(Dat)

% This function identify the range of grid size (partition size, N-by-N)
% based on the inequality satisfying unbiased probability calculation.

% INPUTs

 % Dat: N-by-M matrix for N time recording and M variables 

 % OUTPUTs

% N_min: the smallest grid size

% N_max: the largest grid size satisfying the eq. 17. 

[~, M] = size(Dat);
N_final = zeros(M*(M-1)/2,1);
r=0; 

for k=1:M-1
    for m=k+1:M
     N = 1;
     node1 = Dat(:,k);
     node2 = Dat(:,m);
    
     N_c = 0;
     mean_points = 1;
     
        while (mean_points >= N_c)
        N = N+1;
        [~, ppb] = location([node1, node2],N);        
        
        N_c = sum(ppb(:,2) ~= 0);
        mean_points = length(node1)/N^2;
        end
        r=r+1;
        N_final(r) = N-1; 
    end
end

N_min = round(0.1*mean(N_final))+1;
N_max = round(mean(N_final));
        if N_max <=2 
            N_max = N_max*2;
        end
        
        if N_min <=2 
           N_min = N_min*2;
        end
end