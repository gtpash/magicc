% Forward-Time Centered-Space (FTCS) solver for Fisher-KPP

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


function [A] = assembleA(n, d, h, k)
L = assembleL(n);
L = applyBC(L);
L = (d/h^2)*L;
A = sparse(k*eye(n) + L);
end


function [H] = assembleH(n, k)
H = zeros(n, n^2);
H(:, 1:n+1:end) = eye(n);
H = k*sparse(H);
end


function [N, t] = fwdEuler(A, H, tspan, dt)
% Forward-Euler integrator for "matricized" polynomial system
% NOTE: cell volume fractions are used (theta = 1)
% Inputs:
%   A:          Linera Operator
%   H:          Quadratic Operator
%   tspan:      [t0, tf]
%   dt:         time step
% Outputs:
%   N:          matrix of of solution snapshots
%   t:          solution snapshot times

% simulation times
t = tspan(1):dt:tspan(2);
nt = numel(t);

% matrix to store output
N = zeros(numel(N0), nt);
N(:,1) = N0;

% forward euler time-stepping
for i = 1:nt-1
    N(:, i+1) = N(:, i) + dt*(A*N(:, i) - H*kron(N(:, i), N(:, i)));
end

end


function [N] = proliferation(N, k, theta)
N = k*N.*(1 - N/theta);
end


function [N] = treatment(N, C, alpha1, alpha2, beta1, beta2, t)
    (alpha1*exp(-beta1*t) + alpha2*exp(-beta2*t))*C*N;
end


function [N, t] = FTCS_tx(N0, tspan, txduration, dt, params)
% Forward-Euler integrator for "matricized" polynomial system
% NOTE: cell volume fractions are used (theta = 1)
% Inputs:
%   N0:         initial condition
%   txduration: interval between drug doses
%   tspan:      [t0, tf]
%   dt:         time step
%   params:     structure containing treatment parameters
% Outputs:
%   N:          matrix of of solution snapshots
%   t:          solution snapshot times

% simulation times
t = tspan(1):dt:tspan(2);
nt = numel(t);
nt_tx = numel(0:dt:txduration);

% matrix to store output
N = zeros(numel(N0), nt);
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

function [N, t] = FTCS(L, f, N0, tspan, dt, D, h, k, theta)
% Forward-Euler in Time, Centered Differences in Space
% Solves the governing equation: u_t = D∆u + ku(1-u/theta)
% Inputs:
%   L:          discrete laplacian
%   f:          nonlinear proliferation term [function]
%   tspan:      [t0, tf]
%   k:          proliferation rate
%   theta:      carrying capacity
% Outputs:
%   N:          matrix of of solution snapshots
%   t:          solution snapshot times

% simulation times
t = tspan(1):dt:tspan(2);
nt = numel(t);

% matrix to store output
N = zeros(numel(N0), nt);
N(:,1) = N0;

% forward euler time-stepping
for i = 1:nt-1
    N(:, i+1) = N(:, i) + dt*((D/h^2)*L*N(:, i) + f(N(:, i), k, theta));
end
end
