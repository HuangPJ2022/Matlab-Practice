% this file calculates the steady state values and calls dynare 

% clean up workspace 
clear 
clc

% set parameters 
param_values.beta    = 0.95; 
param_values.delta   = 0.025; 
param_values.alpha   = 0.33;
param_values.gamma     = 4; 
param_values.rho     = 0.8;
param_values.chi     = 1;

% solve for steady state
f = @(x) steadystate_solver(x,param_values); 
steadystate_values = fsolve(f,ones(6,1));  


% define function 
function [y] = steadystate_solver(x,par)  
    % define variables 
    y = zeros(6,1); 
    Z = x(1);
    Y = x(2);
    K = x(3);
    C = x(4);
    I = x(5);
    L = x(6);
    
    % steady state function 
    y(1)       = 1 - Z                                                                 ; 
    y(2)       = K^par.alpha * L^(1-par.alpha) - Y                                     ; 
    y(3)       = (par.alpha/(par.beta^(-1)+par.delta-1))^(1/(1 - par.alpha))*L - K       ; 
    y(4)       = ((1 - par.alpha)/par.chi)*(1 - L)*K^par.alpha*L^(-par.alpha) - C^(par.gamma);           ;
    y(5)       = par.delta*K - I                                                       ; 
    y(6)       = Y - I - C                                                             ;
end 