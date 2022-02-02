%{ 
VirtualDataGeneration: creates virtual tumors with known parameters for building POD operators
Formatting:
    - Struct: named after kp and d values, kp0_000d0_000, decimal replaced with underscore
        - .kp: proliferation used for fwd eval
        - .d: diffusivity used for fwd eval
        - .N_0: initial seeding density
        - .N: Snapshots
        - .t: times of snapshots
    Repeated for each combination of kp and d to be tested
%}

clear all; clc;

%Initialize 50x50 grid
x = 1:1:50;
y = 1:1:50;
[X,Y] = meshgrid(x,y);
x_mid = 25.5; y_mid = 25.5; %Center Point

z = exp(-((X-x_mid).^2+(Y-y_mid).^2)/5);

%Initialize tumor seed with 50% cell density
N_0 = zeros(size(z));
N_0(z>0.05) = 0.5; 
% imagesc(N_0); title('Initial Cell Density'); colorbar; axis image;

clear z x y X Y x_mid y_mid


%Put vector bounds here
t_vec = [2, 4, 6, 8, 10]; %days
kp_vec = [0.1, 0.2]; 
d_vec  = [0.1, 0.5, 1.0];

%Simulate fwd eval for all parameter options
for j = 1:length(d_vec)
    for i = 1:length(kp_vec)
        
        kp = kp_vec(i);
        d  = d_vec(j);
        kp_str = replace(num2str(kp,'%.3f'),'.','_');
        d_str  = replace(num2str(d,'%.3f'),'.','_');
        name = ['kp',kp_str,'d',d_str];
        clear kp_str d_str
        
        eval([name,' = struct;']);
        eval([name,'.kp = kp;']);
        eval([name,'.d  = d;']);
        eval([name,'.t  = t_vec;']);
        eval([name,'.N_0 = N_0;']);
        
        eval([name,'.N = RXDIF_2D(N_0, kp, d, t_vec);']);
        
        %Add noise
%         for k = 1:length(t_vec)
%            eval([name,'.N(:,:,k) = imnoise(',name,'.N(:,:,k),''gaussian'',0,0.001);']); 
%         end
        
        %Testing calibrations
%         eval(['Calibrate_RXDIF(',name,',0.0001, 0.0001);']);
    end
end
clear name i j d kp

%Save workspace to store structs
% save('VirtualData.mat');

