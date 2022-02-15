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


function [H] = assembleH(n)
H = zeros(n, n^2);
H(:, 1:n+1:end) = eye(n);
H = sparse(H);
end


function [N] = proliferation(N, k, theta)
N = k*N.*(1 - N/theta);
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
