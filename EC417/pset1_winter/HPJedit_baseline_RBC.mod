var C Y K L I Z; 
varexo eps; 
parameters beta delta alpha rho gamma chi C_ss Y_ss K_ss L_ss I_ss Z_ss;  

% load steady state and parameter data 
load  dynare_inputs; 

% assign parameter values 
set_param_value('beta'  ,dynare_inputs.param_values.beta  );
set_param_value('delta' ,dynare_inputs.param_values.delta );
set_param_value('alpha' ,dynare_inputs.param_values.alpha );
set_param_value('rho'   ,dynare_inputs.param_values.rho   ); 
set_param_value('gamma'   ,dynare_inputs.param_values.gamma   );
set_param_value('chi'   ,dynare_inputs.param_values.chi   ); 

% load steady state values as parameters
set_param_value('Z_ss'  ,dynare_inputs.ss(1)  );
set_param_value('Y_ss'  ,dynare_inputs.ss(2) );
set_param_value('K_ss'  ,dynare_inputs.ss(3) );
set_param_value('C_ss'  ,dynare_inputs.ss(4)   );
set_param_value('I_ss'  ,dynare_inputs.ss(5)   ); 
set_param_value('L_ss'  ,dynare_inputs.ss(6)   );

model; 
    % variables of interest 
    C^(-gamma)  = beta*(1+alpha*Y(+1)/K-delta)*(C(+1))^(-gamma);
    C^(gamma)  = ((1-alpha)/chi)*(Y/L)*(1-L); 
    Y       = Z*K(-1)^alpha*L^(1-alpha);
    Y       = C + I ;
    K       = (1-delta)*K(-1) + I;
    Z       = Z(-1)^rho*exp(eps); 
end; 

initval; 
    Z = Z_ss;
    Y = Y_ss;
    K = K_ss;
    C = C_ss;
    I = I_ss;
    L = L_ss;
end; 


shocks; 
    var eps = 1^2; 
end; 


stoch_simul(irf=20,order=1,loglinear,irf_plot_threshold=0) C Y K Z I L;

