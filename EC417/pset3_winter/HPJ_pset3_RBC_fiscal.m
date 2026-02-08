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
dynare_inputs.param_values.Phi    = 0.01; 
dynare_inputs.param_values.delta  = 0.025; 
dynare_inputs.param_values.alpha  = 0.33;
dynare_inputs.param_values.sigma  = 0.01; 
dynare_inputs.param_values.chi    = 1;
dynare_inputs.param_values.rho    = 0.8;

% solve for steady state
f_25 = @(x) steadystate_solver(x,dynare_inputs.param_values, 0.25); 
dynare_inputs.ss = fsolve(f_25,ones(10,1));  
f_1 = @(x) steadystate_solver(x,dynare_inputs.param_values, 0.1); 
dynare_inputs.ss_1 = fsolve(f_1,ones(10,1));  
clear f

% save output for dynare 
save('dynare_inputs.mat', 'dynare_inputs'); 

% call dynare file  
dynare HPJedit_baseline_RBC.mod noclearall nolog      

function [y] = steadystate_solver(x,par,tax)  
    y = zeros(10,1); 
    Z = x(1); R = x(2); W = x(3); L = x(4); K = x(5); I = x(6); Y = x(7); C = x(8); G = x(9); T = x(10);

    y(1) = 1 - Z; 
    y(2) = par.Phi + par.delta - R;
    y(3) = (1-par.alpha)*((par.alpha/(par.Phi + par.delta))^(par.alpha/(1 - par.alpha))) - W;
    y(4) = ((1 - par.alpha)/((1 - par.alpha) + (par.chi*(1 - par.delta*(par.alpha/(par.Phi + par.delta)))))) - L;
    y(5) = (par.alpha/(par.Phi + par.delta))^(1/(1 - par.alpha))*L - K; 
    y(6) = par.delta*K - I; 
    y(7) = K^par.alpha * L^(1-par.alpha) - Y; 
    y(8) = ((par.alpha/(par.Phi + par.delta))^(par.alpha/(1-par.alpha)) - par.delta*(par.alpha/(par.Phi + par.delta))^(1/(1-par.alpha)))*L - C*(1 + T);
    y(9) = T*C - G;
    y(10) = tax - T;
end