% this file calculates the steady state values and calls dynare 

% for this file to work, you must 
% (1) install dynare 
% (2) add the dynare folder to the matlab paths 
% (3) save this file, as well as baseline_RBC.mod, baseline_RBC_linear.mod
%     to the same folder. The file produces impulse response functions
%     using first order perturbation in Dynare as discussed in the lecture.

clear 
clc

% set parameters 
dynare_inputs.param_values.beta    = 0.95; 
dynare_inputs.param_values.delta   = 0.025; 
dynare_inputs.param_values.alpha   = 0.33;
dynare_inputs.param_values.gamma     = 4; 
dynare_inputs.param_values.chi     = 1;

rho_values = 0:0.2:1.2;

% 2. Loop through them
for i = 1:length(rho_values)
    
    % Assign the current value to your structure
    fprintf('Running iteration %d with rho = %.1f\n', i, rho_values(i));
    dynare_inputs.param_values.rho = rho_values(i);
    
    % solve for steady state
    f = @(x) steadystate_solver(x,dynare_inputs.param_values); 
    dynare_inputs.ss = fsolve(f,ones(6,1));  
    clear f
    
    % save output for dynare 
    save('dynare_inputs.mat', 'dynare_inputs'); 
    
    % call dynare file  
    dynare baseline_RBC.mod noclearall nolog         %levels file
    %dynare baseline_RBC_linear.mod                  %linear file 

end

function [y] = steadystate_solver(x,par)  
    y = zeros(6,1); 
    Z = x(1); Y = x(2); K = x(3); C = x(4); I = x(5); L = x(6);
    
    y(1) = 1 - Z; 
    y(2) = K^par.alpha * L^(1-par.alpha) - Y; 
    y(3) = (par.alpha/(par.beta^(-1)+par.delta-1))^(1/(1 - par.alpha))*L - K; 
    y(4) = ((1 - par.alpha)/par.chi)*(1 - L)*K^par.alpha*L^(-par.alpha) - C^(par.gamma);
    y(5) = par.delta*K - I; 
    y(6) = Y - I - C;
end