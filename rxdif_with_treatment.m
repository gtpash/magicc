function [L] = assembleL(n)
I = eye(n);
e = ones(n,1);
T = spdiags([e -4*e e],[-1 0 1],n,n);
S = spdiags([e e],[-1 1],n,n);
L = (kron(I,T) + kron(S,I));
end


function [L] = applyBC(L)
% Apply Neumann BC discrete laplacian

n = sqrt(size(L,1));

% Apply BC to left wall
for i = 1:n:n^2
    L(i, i+1) = 2;
end

% Apply BC to right wall
for i = 1:n:n^2
    L(i, i+1) = 2;
end

% Apply BC to bottom wall
for i = 1:n
    L(i, i+n) = 2;
end

% Apply BC to top wall
for i = 1:n
    L(end-i, end-i-n) = 2;
end
end


function [H] = assembleH(n, k)
H = zeros(n, n^2);
H(:, 1:n+1:end) = eye(n);
H = k*sparse(H);
end


function [N] = treatment(N, C, alpha1, alpha2, beta1, beta2, t)
    (alpha1*exp(-beta1*t) + alpha2*exp(-beta2*t))*C*N;
end


function [N, t] = FTCS_tx(N0, params)
% Forward-Euler integrator for "matricized" polynomial system
% NOTE: cell volume fractions are used (theta = 1)
% Inputs:
%   N0:         initial condition
%   params:     structure containing treatment parameters
% Outputs:
%   N:          matrix of of solution snapshots
%   t:          solution snapshot times

% unpack parameters
tspan = params.tspan;
dt = params.dt;
txduration = params.txduration;
alpha1 = params.alpha1;
alpha2 = params.alpha2;
beta1 = params.beta1;
beta2 = params.beta2;

% simulation times
t = tspan(1):dt:tspan(2);
nt = numel(t);
nt_tx = numel(0:dt:txduration);

% matrix to store output
n = numel(N0);
N = zeros(n, nt);
N(:,1) = N0;

% form operators
L = assembleL(n);
L = applyBC(L);
H = assembleH(n, k);

% forward euler time-stepping
for i = 1:nt-1
    N(:, i+1) = N(:, i) + dt*(d/h^2)*L*N(:, i) + dt*k*eye(n)*N(:, i) ...
        - dt*H*kron(N(:, i), N(:, i)) ...
        - dt*treatment(N, C, alpha1, alpha2, beta1, beta2, mod(nt, nt_tx));
end
end
