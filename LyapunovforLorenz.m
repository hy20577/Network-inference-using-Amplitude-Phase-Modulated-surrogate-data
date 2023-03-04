%% Lyapunov Exponent for Lorenz Equations with the following arguments:
%[ tvalues, Lyapunov_spectrum, Orbits] = LyapunovforLorenz(x0, t0, tf, k, coupling_strength, sigma, rho, beta, Adj)
% x0 --> initial conditions
% t0 --> initial time
% tf --> integration time
% k --> number of iteration
% coupling_strength --> coupling strength 
% sigma, rho, beta --> Parameters of Lorenz (as sigma =10, rho=28, beta=8/3)
% Adj --> Adjacency Matrix of the initial Network

function [tvalues, Lyapunov_spectrum, Orbits] = LyapunovforLorenz(x0, t0, tf, k, coupling_strength, sigma, rho, beta, Adj)

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
        [T, Y] = ode45(@(t,y) lorenz_ext(t,y, Adj, coupling_strength, sigma, rho, beta), tspan, y0, options);
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

function f=lorenz_ext(t,X, Adj, K, sigma, rho, beta)

M = size(Adj,2);
% Parameters
% sigma = 16;  %10;
% rho = 45.92;
% beta = 4; %8/3;

x=X(1:M); y=X(M+1:2*M); z=X(2*M+1:3*M);

f=zeros(3*M+9*M^2,1);

        for i=1:M
            sum1 = 0;
                for j=1:M
                    sum1 = sum1 + Adj(i,j)*(x(j)-x(i));
                end
            f(i) = sigma*(y(i) - x(i)) + K*sum1;
            f(M+i) = x(i)*(rho-z(i))-y(i);
            f(2*M+i) = x(i)*y(i)-beta*z(i);
        end

x_dx = zeros(M,M);
    for i=1:M
        for j=1:M
            if i ==j
            x_dx(i,j) = -sigma-K*sum(Adj(i,:));
            else
                x_dx(i,j) = K*Adj(i,j);
            end
        end
    end

 x_dy = zeros(M,M);
     for i=1:M
         x_dy(i,i)=sigma;
     end
 x_dz = zeros(M,M);
 y_dx = zeros(M,M);
     for i=1:M
             y_dx(i,i) = rho-z(i);
     end
y_dy = zeros(M,M);
    for i=1:M
    y_dy(i,i)= -1;
    end
    y_dz = zeros(M,M);
    for i=1:M
    y_dz(i,i)= -x(i);
    end
z_dx = zeros(M,M);

    for i=1:M
    z_dx(i,i)= y(i);
    end

z_dy = zeros(M,M);

    for i=1:M
    z_dy(i,i)= x(i);
    end
z_dz = zeros(M,M);

    for i=1:M
    z_dz(i,i)= -beta;
    end

    Jacob_Lorenz = [x_dx, x_dy, x_dz; y_dx, y_dy, y_dz; z_dx, z_dy, z_dz];

Y = zeros(3*M,3*M);
for i=1:3*M
Y(:,i) = X(i*3*M+1:(i+1)*3*M);
end

f(3*M+1:3*M+(3*M)^2)= Jacob_Lorenz*Y;

end


