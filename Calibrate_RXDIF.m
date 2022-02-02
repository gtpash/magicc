%{ 
Calibrate_RXDIF: Levenberg-Marquardt fitting for the 2D reaction diffusion equation
Slow! Builds jacobian for every change of parameters

Inputs:
    - tumor: struct containing all tumor parameters
        - N_0: initial cell count
        - N: snapshots of growth
        - t: time of snapshots
        - kp: true proliferation
        - d: true diffusivity
    - kp_g: proliferation guess
    - d_g: diffusivity guess
Outputs:
    - kp_cal: calibrated proliferation estimate
    - d_cal: calibrated diffusivity estimate
    - MSE: data mismatch at convergence or max iterations

%}

function [kp_cal, d_cal, MSE] = Calibrate_RXDIF(tumor, kp_g, d_g)
    start = tic;
    
    %LM parameters
    e_tol    = 1e-5;
    e_conv   = 1e-8;
    max_it   = 500;
    delta    = 1.001;
    pass     = 7;
    fail     = 9;
    alpha    = 1;
    j_change = 1;
    
    %Setup parameter bounds and guesses
    kp_curr = kp_g;
    d_curr  = d_g;
    kp_up  = 10;
    kp_low = 1e-6;
    d_up   = 10;
    d_low  = 1e-6;
    
    
    %Pull out struct variables
    N_0 = tumor.N_0;
    N = tumor.N;
    t = tumor.t;
    [sy,sx] = size(N_0);
    mid = round(sy/2);
    
    %initialize SSE
    N_sim = RXDIF_2D(N_0, kp_curr, d_curr, t);
    SSE   = sum((N - N_sim).^2,'all');
    
    %Calibration Loop
    iteration = 1;
    while(iteration<=max_it && SSE>e_tol)
        %Build Jacobian
        if(j_change==1)
            %Perturb and find difference
            kp_p = kp_curr*delta;
            d_p  = d_curr*delta;
            dif_kp = kp_p - kp_curr;
            dif_d  = d_p - d_curr;
            
            %Simulate perturbed params
            N_kpp = RXDIF_2D(N_0, kp_p, d_curr, t);
            N_dp  = RXDIF_2D(N_0, kp_curr, d_p, t);
            
            %Fill J
            J = zeros(numel(N), 2);
            J(:,1) = reshape(N_kpp - N_sim, [], 1)/dif_kp;
            J(:,2) = reshape(N_dp - N_sim, [], 1)/dif_d;
        end
        
        %Calculate update with current regularization
        residuals = reshape(N - N_sim, [], 1);
        [update,flags] = bicgstab((J'*J + alpha*diag(diag(J'*J))),(J'*residuals),1e-5,250);
        
        %Update parameters to test
        kp_test = kp_curr + update(1);
        if(kp_test<kp_low)
            kp_test = kp_low;
        elseif(kp_test>kp_up)
            kp_test = kp_up;
        end
        d_test  = d_curr + update(2);
        if(d_test<d_low)
            d_test = d_low;
        elseif(d_test>d_up)
            d_test = d_up;
        end
        
        %Test SSE calculation
        N_test = RXDIF_2D(N_0, kp_test, d_test, t);
        SSE_test = sum((N - N_test).^2,'all');
        
        %Evaluate new error
        if(SSE_test<SSE)
            N_sim = N_test;
            kp_curr = kp_test;
            d_curr  = d_test;
            
            if(SSE-SSE_test < e_conv)
                disp(['Algorithm converged on iteration: ',num2str(iteration)]);
                break;
            end
            
            SSE = SSE_test;
            
            alpha = alpha/pass;
            if(SSE<e_tol)
                disp(['Tolerance reached on iteration: ',num2str(iteration)]);
                break;
            end
        else
            alpha = alpha*fail;
        end
        iteration = iteration + 1;
    end
    stop = toc(start);
    disp(['Time to fit: ',num2str(stop/60),' min']);
    
    kp_cal = kp_curr;
    d_cal  = d_curr;
    MSE    = mean((N-N_sim).^2,'all');
    
    
    %Visualize and print fits, comment out later
%     colors = ['r','b','g','y','m'];
%     x = 1:sx;
%     figure
%     plot(N_0(mid,:),'k-','DisplayName','Initial','LineWidth',2);
%     hold on
%     for i = 1:length(t)
%         c = colors(1); colors(1) = [];
%         meas_str = ['Meas Day ',num2str(t(i))];
%         sim_str = ['Fit Day ',num2str(t(i))];
%         scatter(x,squeeze(N(mid,:,i)), 10, c, 'DisplayName', meas_str);
%         plot(squeeze(N_sim(mid,:,i)),'-','DisplayName',sim_str,'Color',c, 'LineWidth', 2);
%     end
%     title(['Kp = ',num2str(tumor.kp),'; D = ',num2str(tumor.d)]);
%     legend();
%     
%     kp_err = 100*(kp_cal - tumor.kp)/tumor.kp;
%     d_err  = 100*(d_cal - tumor.d)/tumor.d;
%     disp(['Kp error: ', num2str(kp_err)]);
%     disp(['D error: ', num2str(d_err)]);
    
end