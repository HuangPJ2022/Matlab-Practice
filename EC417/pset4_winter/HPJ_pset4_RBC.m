clear 
clc

% set parameters 
dynare_inputs.param_values.beta   = 0.95;
dynare_inputs.param_values.delta0 = 0.025; 
dynare_inputs.param_values.alpha  = 0.33;
dynare_inputs.param_values.sigma  = 0.01; 
dynare_inputs.param_values.omega  = ((1/0.95) - 1 + 0.025) / 0.025;
dynare_inputs.param_values.rho    = 0.8;
dynare_inputs.param_values.gamma  = 1;

% solve for steady state
f = @(x) baseline_steadystate_solver(x,dynare_inputs.param_values); 
dynare_inputs.ss_baseline = fsolve(f,ones(10,1));
clear f

f = @(x) compare_steadystate_solver(x,dynare_inputs.param_values); 
dynare_inputs.ss_compare = fsolve(f,ones(10,1));
clear f

% save output for dynare 
save('dynare_inputs.mat', 'dynare_inputs'); 

% call dynare file  
dynare HPJedit_baseline_RBC.mod noclearall nolog 
dynare HPJedit_compare_RBC.mod noclearall nolog

function [y] = baseline_steadystate_solver(x,par)  
    y = zeros(10,1); 
    Z = x(1); R = x(2); W = x(3); L = x(4); K = x(5); I = x(6); Y = x(7); C = x(8); V = x(9); d = x(10);

    y(1) = 1 - Z; 
    y(2) = d + (1/par.beta) -1 - R;
    y(3) = (1-par.alpha)*((par.alpha/R)^(par.alpha/(1 - par.alpha))) - W;
    y(4) = (1-par.alpha)/(1-(d*par.alpha/R)) - L^(1+(1/par.gamma));
    y(5) = (par.alpha/R)^(1/(1 - par.alpha))*L - K; 
    y(6) = d*K - I; 
    y(7) = K^par.alpha * L^(1-par.alpha) - Y; 
    y(8) = (1-(d*par.alpha)/R)*(par.alpha/R)^(par.alpha/(1-par.alpha))*L - C;
    y(9) = 1 - V;
    y(10) = par.delta0 - d;
end

function [y] = compare_steadystate_solver(x,par)  
    y = zeros(10,1); 
    Z = x(1); R = x(2); W = x(3); L = x(4); K = x(5); I = x(6); Y = x(7); C = x(8); V = x(9); d = x(10);

    y(1) = 1 - Z; 
    y(2) = d + (1/par.beta) -1 - R;
    y(3) = (1-par.alpha)*((par.alpha/R)^(par.alpha/(1 - par.alpha))) - W;
    y(4) = W^(par.gamma) - L;
    y(5) = (par.alpha/R)^(1/(1 - par.alpha))*L - K; 
    y(6) = d*K - I; 
    y(7) = K^par.alpha * L^(1-par.alpha) - Y; 
    y(8) = (1-(d*par.alpha)/R)*(par.alpha/R)^(par.alpha/(1-par.alpha))*L - C;
    y(9) = 1 - V;
    y(10) = par.delta0 - d;
end