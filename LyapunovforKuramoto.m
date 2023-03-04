%% Lyapunov Exponent for Deterministic Kuramoto System with the following arguments:

% x0 --> initial conditions
% t0 --> initial time
% tf --> integration time
% k --> number of iteration
% coupling_strength --> coupling strength 
% Adj --> Adjacency Matrix of the initial Network
% w --> M-by-1 array corresponding the internal frequency of each
% oscillator.
function [tvalues, Lyapunov_spectrum, Orbits] = LyapunovforKuramoto(x0, t0, tf, k, coupling_strength, Adj, w)

%%% ode45 function is used to solve ode equation, needed to revise for
%%% other solution.

% ext_odefun - extended ode system needed as input


% x0: initial condition for the ode system
% t0: initial time
% tf: final time
% k: number of iteration
% How many data you want to discard in computation of LE as 

n = length(x0); % dimension of extended ode function.
stept = (tf-t0)/k;
sum1= 0;

% Setting initial conditions
t= t0;
y0 = zeros(1, n^2+n);
y0(1:n) = x0;

for i=1:n
    y0((n+1)*i) = 1; % initial deviation vectors: unit vectors.
end

Orbits = zeros(k,n);
tvalues = zeros(k,1);
evolution_LYA = zeros(k,n+1);

    for ITERLAP=1:k   %Main loop to compute Lyapunov
        
        tspan= [t t+stept];
        options = odeset('RelTol',1e-8);
        [T, Y] = ode45(@(t,y) Kuramoto_extended(t, y, coupling_strength, Adj, w), tspan, y0, options);
        tvalues(ITERLAP) = t+stept;
        Y = Y(end,:);
        y0(1:n) = Y(1:n);
        Orbits(ITERLAP,:) = Y(1:n);
        
        K = reshape(Y(n+1:n^2+n),n,n);
        
%         for j=1:n
%      
%             K(:,j) = Y(n*j+1 : n*(j+1));  %% create deviation matrix
%             
%         end
        
        K = gramschmidt(K);    %% Orthogonalization
        
        norm_K = zeros(1, n);
        NK = zeros(n,n);
        
        for j=1:n
            norm_K(j) = norm(K(:,j));
            NK(:,j) = K(:,j)/norm_K(j);
        end
       
        sum1 = sum1 + log(norm_K);
       

       evolution_LYA(ITERLAP,1) = t; 
       evolution_LYA(ITERLAP,2:end) = sum1/(t-t0);

        t = t+stept;
        
        y0(n+1:n^2+n) = reshape(NK,1,n^2);
        
%         for j=1:n
%          y0(n*j+1 : n*(j+1)) = K(:,j)/norm_K(j);  %% prepare deviation matrix for next iteration.   
%         end
   
    if mod(ITERLAP,1e2) == 0
        fprintf('Progression: %s %% \n', num2str(ITERLAP*100/k))
    end
    
    end
Lyapunov_spectrum = sum1/(tf-t0);
Lyapunov_spectrum = sort(Lyapunov_spectrum, 'descend');
% writematrix(evolution_LYA, 'Time_evolution_Lyapunov_Spec.txt')

end

function U = gramschmidt(V)
% GRAM_SCHMIDT - Classic Gram-Schmidt Process
%   Input vectors should be the columns of input matrix.
%   Output = unitary, orthogonal vectors in columns of matrix

n = size(V, 2); % number of columns

U(:,1) = V(:,1);

% find next orthogonal column
for ii = 2:n
    b = V(:,ii);
    c = b;
    for k = 1:ii-1
        a = U(:,k);
        c = c - (a*a'*b)/(a'*a);
    end
    U(:,ii) = c;
end
end

%% Functions 


function theta = Kuramoto_extended(t, y, K, A, w)

% A adjacency Matrix
% K is the coupling strength of the model
% y is extended inputs. First M is the initial conditions for M
% oscillators, next ones deviation vectors. 

M = size(A,2); % number of oscillators
theta = zeros(M+M^2,1);
Jacob_mat = zeros(M,M); 
% Computation of orbits and Jacobian

    for i=1:M
        sum1 = 0 ;
        sum2 = 0;
        for j=1:M
            sum1 = sum1+ A(i,j)*sin(y(j)-y(i));
            sum2 = sum2 + A(i,j)*cos(y(j)-y(i));
            if i ~= j
                Jacob_mat(i,j) = (K/M)*A(i,j)*cos(y(j)-y(i));
            end

        end
    theta(i) = w(i) + (K/M)*sum1;  
    Jacob_mat(i,i) = -(K/M)*sum2;
    end


Y = zeros(M,M);
for i=1:M
Y(:,i) = y(i*M+1:(i+1)*M);
end

theta(M+1:M+M^2)= Jacob_mat*Y;


end