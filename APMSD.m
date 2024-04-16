function SD_overall = APMSD(Dat, k, varargin)
   
    if nargin<3
        pc1 = 1;
        pc2 = 0.5;
    elseif nargin == 3 
        pc1 = varargin{1}; % pc1; 1st argument of the varargin.
        pc2 = 0.5;
    elseif nargin == 4
        pc1 = varargin{1}; % pc1; 1st argument of the varargin.
        pc2 = varargin{2}; % pc2; 2nd argument of the varargin.
    end

        M = size(Dat,2);
        n = size(Dat,1);
        
        AS = hilbert(Dat); % analytic signal
        A = abs(AS); 
        ph = angle(AS); % phases
        
        SD_overall = zeros(n,M,k);
        mxph = max(max(ph));
        minph = min(min(ph));
        
        mxA = max(max(A));
        mnA = min(min(A));


        for i=1:k
            perm1 = randperm(size(A,1));
            indices_A = perm1(1:ceil(pc1*size(A,1)));
            A_sd = A; 
            A_sd(indices_A, :) = rand(length(indices_A), M) * (mxA-mnA)+mnA; 

            perm2 = randperm(size(A,1));
            indices_ph = perm2(1:ceil(pc2*size(A,1)));

            theta_sd = ph;
            theta_sd(indices_ph, :) = rand(length(indices_ph),M)* (mxph- minph)+minph;
            SD_overall(:,:,i) = real(A_sd.*exp(1i*theta_sd));
        end

    
end







