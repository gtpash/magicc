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

%% Generate Synthetic Tumor Data (no treatment)

clear; clc; close all;

%  parameter bounds
kp_vec = linspace(0.001, 0.1, 10);
d_vec = linspace(1e-6, 1e-3, 10);
[kps, ds] = meshgrid(kp_vec, d_vec);
kd = [kps(:), ds(:)];  % pairs of (k, d)
nmodes = 0*kps;  % store how many required modes

CNRG_THRESHOLD = 0.99;

n = 50;

viewN = @(N, nmesh) imagesc(reshape(N, nmesh, nmesh));

% initializations
N0 = initialize_tumor(n);
params.dt = 0.1;    % time step [day]
params.h = 1;       % spatial discretization [mm]
params.tspan = [0, 56];

% nsnap = params.tspan(end) / params.dt;
snapmat = [];

for ii = 1:length(kp_vec)
    for jj = 1:length(d_vec)
        params.k = kp_vec(ii);
        params.d = d_vec(jj);

        [N, t] = FTCS_tx(N0, params);
        snapmat = [snapmat; N];

        [sv, cnrg] = computesv(N);

        n_cnrg = numel(cnrg(cnrg < CNRG_THRESHOLD));
        nmodes(ii, jj) = n_cnrg + 1;
    end
end

% fprintf("Modes required for Cumulative Energy to be %.3f:\t%i\n", CNRG_THRESHOLD, n_cnrg+1);

%% Plot Singular Value Decay + Cumulative Energy

[sv, cnrg] = computesv(snapmat);
n_cnrg = numel(cnrg(cnrg < CNRG_THRESHOLD));

figure;
semilogy(sv/max(sv), '-ok', 'linewidth', 2);
hold on; grid on;
semilogy(sv(1:n_cnrg)/max(sv), 'or', 'linewidth', 2);
xlim([-10, 600]);
rectangle('Position', [-5, .01, 20, 50], 'linewidth', 2', 'linestyle', '--');
title('Singular Value Decay of Snapshot Matrix', 'fontsize', 22);
xlabel('Mode, i', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Singular Value, $\sigma_i/\|\sigma\|_\infty$', 'interpreter', 'latex', 'fontsize', 20);

figure;
semilogy(sv(1:15)/max(sv),'-ok','linewidth',2)
hold on; grid on;
semilogy(sv(1:n_cnrg)/max(sv), 'or', 'linewidth', 2);
axis([-5, 20, .01, 50]);
ax = gca;
ax.LineWidth = 2;

figure;
plot(cnrg,'-ok','linewidth',2);
hold on; grid on;
plot(cnrg(1:n_cnrg),'or','linewidth',2)
title('Cumulative Energy', 'fontsize', 22)
ylabel('Cumulative Energy', 'FontSize', 20);
xlabel('Mode, i', 'FontSize', 20);
rectangle('Position', [-5, 0.8, 50, 0.2], 'linewidth', 2', 'linestyle', '--');
xlim([-10, 600]);
ylim([0, 1]);

figure;
plot(cnrg,'-ok','linewidth',3)
hold on; grid on;
plot(cnrg(1:n_cnrg),'or','linewidth',3)
axis([-5, 50 0.8, 1]);
ax = gca;
ax.LineWidth = 2;

%% Generate Tumor Data with ChemoTherapy Treatment Terms
clear; clc; close all;

viewN = @(N, nmesh) imagesc(reshape(N, nmesh, nmesh));

PATIENT_DATA_FILE = "../data/pcr/middle_slice_data_patient290.mat";
load(PATIENT_DATA_FILE);

CNRG_THRESHOLD = 0.99;

%  parameter bounds
kp_vec = linspace(0.001, 0.1, 10);
d_vec = linspace(1e-6, 1e-3, 10);
[kps, ds] = meshgrid(kp_vec, d_vec);
kd = [kps(:), ds(:)];  % pairs of (k, d)
nmodes = 0*kps;  % store how many required modes

% initial condition
n = 26;             % number of nodes in each dimension
% N0 = NTC_Visit1;
% N0 = N0(:);
N0 = initialize_tumor(n);

params.dt = 0.1;    % time step [day]
params.h = 1;       % spatial discretization [mm]
params.tspan = [0, 56];
params.use_tx = true;

params.txduration = 56/2;
params.alpha1 = 0.5;
params.alpha2 = 0.8;
params.beta1 = 0.3;
params.beta2 = 1;
params.C = AUC(:);

snapmat = [];

for ii = 1:length(kp_vec)
    for jj = 1:length(d_vec)
        params.k = kp_vec(ii);
        params.d = d_vec(jj);

        [N, t] = FTCS_tx(N0, params);
        snapmat = [snapmat; N];

        [sv, cnrg] = computesv(N);

        n_cnrg = numel(cnrg(cnrg < CNRG_THRESHOLD));
        nmodes(ii, jj) = n_cnrg + 1;
    end
end

[sv, cnrg] = computesv(snapmat);
plotsv(sv, cnrg);

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
for i = n:n:n^2
    L(i, i-1) = 2;
end

% Apply BC to bottom wall
for i = 1:n
    L(i, i+n) = 2;
end

% Apply BC to top wall
for i = 0:n
    L(end-i, end-i-n) = 2;
end
end


function [H] = assembleH(n, k)
H = sparse(zeros(n, n^2));
H(:, 1:n+1:end) = eye(n);
H = k*sparse(H);
end


function [N] = apply_treatment(N, C, alpha1, alpha2, beta1, beta2, t)
    N = (alpha1*exp(-beta1*t) + alpha2*exp(-beta2*t))*C.*N;
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
if isfield(params, 'use_tx')
    use_tx = params.use_tx;
else
    use_tx = false;
end
if use_tx
    txduration = params.txduration;
    alpha1 = params.alpha1;
    alpha2 = params.alpha2;
    beta1 = params.beta1;
    beta2 = params.beta2;
    nt_tx = numel(0:dt:txduration);
    C = params.C;
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
    if use_tx
        treat = dt*apply_treatment(N(:,i), C, alpha1, alpha2, beta1, beta2, mod(i, nt_tx));
    else
        treat = 0;
    end

    N(:, i+1) = N(:, i) + dt*(d/h^2)*L*N(:, i) + proliferation(N(:, i), k, theta) - treat;
%     N(:, i+1) = N(:, i) + dt*(d/h^2)*L*N(:, i) + dt*k*eye(n)*N(:, i) ...
%         - dt*H*kron(N(:, i), N(:, i)) - treat;

end
end


function [sv, cnrg] = computesv(N)
    sv = svd(N);

    cnrg = cumsum(sv.^2);
    cnrg = cnrg/cnrg(end);
end

function [cnrg] = plotsv(sv, cnrg)
    figure; 
    semilogy(sv/max(sv), 'k', 'linewidth', 3);
    title('Singular Value Decay of Snapshot Matrix', 'fontsize', 22);
    xlabel('i', 'interpreter', 'latex', 'fontsize', 16);
    ylabel('$\sigma_i/\|\sigma\|_\infty$', 'interpreter', 'latex', 'fontsize', 16);

    figure;
    plot(cnrg, 'kx', 'LineWidth', 3);
    title('Cumulative Energy', 'fontsize', 22)
    ylabel('Cumulative Energy', 'interpreter', 'latex', 'FontSize', 16);
    ylim([0, 1]);
    xlabel('$i$', 'interpreter', 'latex', 'FontSize', 16);
end

