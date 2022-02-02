%{ 
RXDIF_2D: Forward evaluation of reaction diffusion equation in 2D
Inputs:
    - initial: initial cell density map
    - kp: global proliferation rate (1/day)
    - d: global diffusivity (mm^2/day)
    - t: time array to output solution vector (days)
Outputs:
    - N_sim: cell densitity maps at times specified by t
%}

function [N_sim] = RXDIF_2D(initial, kp, d, t)
    %FDM grid spacing
    dx = 1;
    dy = 1;
    dt = 0.1;
    
    theta = 1; %If using volume fractions

    t_ = (t/dt) + 1; %Indices of densities to output
    
    %Intialize solution matrix
    [sy,sx] = size(initial);
    st = length([0:dt:t(end)]);
    Sim = zeros(sy,sx,st);
    Sim(:,:,1) = initial;
    
    %Time stepping
    for k = 1:st-1
        temp = zeros(sy,sx);
        N = Sim(:,:,k);
        
        %Space stepping
        for i = 1:(sy*sx)
            
            %Boundary Check
            [y,x] = ind2sub([sy,sx],i);
            boundary = [0,0]; %[y x z], exists on -1 to 1
            try %negative y-check
                test = temp(y-1,x);
            catch
                boundary(1) = -1; 
            end
            try %positive y-check
                test = temp(y+1,x);
            catch
                boundary(1) = 1; 
            end
            
            try %negative x-check
                test = temp(y,x-1);
            catch
                boundary(2) = -1; 
            end
            try %positive x-check
                test = temp(y,x+1);
            catch
                boundary(2) = 1; 
            end
            clear test
            
            %FDM in Y direction
            if(boundary(1)==0)
                inv_y = d*(N(i+1)-2*N(i)+N(i-1))/(dy^2);
                
            elseif(boundary(1)==1)
                inv_y = d*(-2*N(i)+2*N(i-1))/(dy^2);
                
            elseif(boundary(1)==-1)
                inv_y = d*(-2*N(i)+2*N(i+1))/(dy^2);
                
            end
            
            %FDM in X direction
            if(boundary(2)==0)
                inv_x = d*(N(i+sy)-2*N(i)+N(i-sy))/(dx^2);
                
            elseif(boundary(2)==1)
                inv_x = d*(-2*N(i)+2*N(i-sy))/(dx^2);
                
            elseif(boundary(2)==-1)
                inv_x = d*(-2*N(i)+2*N(i+sy))/(dx^2);
                
            end
            
            invasion = inv_y + inv_x;
            prolif   = N(i)*kp*(1-(N(i)/theta));
            
            temp(i) = N(i) + dt*(invasion + prolif);
        end
        
        Sim(:,:,k+1) = temp;
    end
    
    N_sim = Sim(:,:,t_);
end