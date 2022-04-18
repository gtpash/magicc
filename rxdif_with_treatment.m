%{ 
rxdif_with_treatment: creates virtual tumors with known parameters for building POD operators
Formatting:
    - Struct: named after kp and d values, kp0_000d0_000, decimal replaced with underscore
        - .kp: proliferation used for fwd eval
        - .d: diffusivity used for fwd eval
        - .N_0: initial seeding density
        - .N: Snapshots
        - .t: times of snapshots
    Repeated for each combination of kp and d to be tested
%}

%% Initializations

clear; clc; close all;

%  parameter bounds
kp_vec = linspace(0.001, 0.1, 10);
d_vec = linspace(1e-6, 1e-3, 100);

n = 50;

viewN = @(N) imagesc(reshape(N, n, n));

% initializations
N0 = initialize_tumor(n);
params.dt = 0.1;    % time step [day]
params.h = 1;       % spatial discretization [mm]
params.tspan = [0, 56];
params.k = kp_vec(2);
params.d = d_vec(5);

[N, t] = FTCS_tx(N0, params);

cnrg = plotsv(N);
CNRG_THRESHOLD = 0.99;
n_cnrg = numel(cnrg(cnrg < CNRG_THRESHOLD));
fprintf("Modes required for Cumulative Energy to be %.3f:\t%i\n", CNRG_THRESHOLD, n_cnrg+1);

%% Generate Tumor Data


% loop through (d, kp) combos
% for jj = 1:length(d_vec)
%     for ii = 1:length(kp_vec)
%         kp = kp_vec(ii);
%         d = d_vec(jj);
% 
%         kp_str = replace(num2str(kp,'%.3f'),'.','_');
%         d_str  = replace(num2str(d,'%.3f'),'.','_');
%         name = ['kp',kp_str,'d',d_str];
% 
%         
% 
%     end
% end


%% Helper Function Definitions

function [N0] = initialize_tumor(n)
%Initialize 50x50 grid
x = 1:1:n;
y = 1:1:n;
[X,Y] = meshgrid(x,y);
x_mid = (n+1)/2; y_mid = x_mid; %Center Point

z = exp(-((X-x_mid).^2+(Y-y_mid).^2)/5);

%Initialize tumor seed with 50% cell density
N0 = zeros(size(z));
N0(z>0.05) = 0.5;
N0 = N0(:);
end

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
H = sparse(zeros(n, n^2));
H(:, 1:n+1:end) = eye(n);
H = k*sparse(H);
end


function [N] = treatment(N, C, alpha1, alpha2, beta1, beta2, t)
    (alpha1*exp(-beta1*t) + alpha2*exp(-beta2*t))*C*N;
end


function [N] = proliferation(N, k, theta)
N = k*N.*(1 - N/theta);
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
h = params.h;
k = params.k;
d = params.d;
theta = 1;
if isfield(params, 'treatment')
    treatment = params.treatment;
else
    treatment = false;
end
if treatment
    txduration = params.txduration;
    alpha1 = params.alpha1;
    alpha2 = params.alpha2;
    beta1 = params.beta1;
    beta2 = params.beta2;
    nt_tx = numel(0:dt:txduration);
end

% simulation times
t = tspan(1):dt:tspan(2);
nt = numel(t);

% matrix to store output
N = zeros(numel(N0), nt);
N(:,1) = N0;

% form operators
L = assembleL(sqrt(numel(N0)));
L = applyBC(L);
% H = assembleH(n, k);

% forward euler time-stepping
for i = 1:nt-1
    if treatment
        treat = dt*treatment(N, C, alpha1, alpha2, beta1, beta2, mod(nt, nt_tx));
    else
        treat = 0;
    end

    N(:, i+1) = N(:, i) + dt*(d/h^2)*L*N(:, i) + proliferation(N(:, i), k, theta) - treat;
%     N(:, i+1) = N(:, i) + dt*(d/h^2)*L*N(:, i) + dt*k*eye(n)*N(:, i) ...
%         - dt*H*kron(N(:, i), N(:, i)) - treat;    

end
end

function [cnrg] = plotsv(N)
    sv = svd(N);
    figure; 
    semilogy(sv/max(sv), 'k', 'linewidth', 3);
    title('Singular Value Decay of Snapshot Matrix');
    xlabel('i');
    ylabel('$\sigma_i/\|\sigma\|_\infty$', 'interpreter', 'latex');

    figure;
    cnrg = cumsum(sv.^2);
    cnrg = cnrg/cnrg(end);
    plot(cnrg, 'kx', 'LineWidth', 3);
    ylabel('Cumulative Energy', 'interpreter', 'latex', 'FontSize', 16);
    ylim([0, 1]);
    xlabel('$i$', 'interpreter', 'latex', 'FontSize', 16);
end
