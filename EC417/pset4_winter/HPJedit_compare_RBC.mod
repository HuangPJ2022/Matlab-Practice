var Z R W L K I Y C V d ;
varexo eps; 
parameters beta delta0 sigma alpha rho omega gamma Z_ss R_ss W_ss L_ss K_ss I_ss Y_ss C_ss V_ss d_ss;

% load steady state and parameter data 
load  dynare_inputs; 

% assign parameter values 
set_param_value('beta'  ,dynare_inputs.param_values.beta  );
set_param_value('delta0' ,dynare_inputs.param_values.delta0 );
set_param_value('sigma' ,dynare_inputs.param_values.sigma );
set_param_value('alpha'   ,dynare_inputs.param_values.alpha   ); 
set_param_value('rho'   ,dynare_inputs.param_values.rho   );
set_param_value('omega'   ,dynare_inputs.param_values.omega   ); 
set_param_value('gamma'   ,dynare_inputs.param_values.gamma   ); 

% load steady state values as parameters
set_param_value('Z_ss'  ,dynare_inputs.ss_compare(1)  );
set_param_value('R_ss'  ,dynare_inputs.ss_compare(2) );
set_param_value('W_ss'  ,dynare_inputs.ss_compare(3) );
set_param_value('L_ss'  ,dynare_inputs.ss_compare(4)   );
set_param_value('K_ss'  ,dynare_inputs.ss_compare(5)   ); 
set_param_value('I_ss'  ,dynare_inputs.ss_compare(6)   );
set_param_value('Y_ss'  ,dynare_inputs.ss_compare(7)   );
set_param_value('C_ss'  ,dynare_inputs.ss_compare(8)   );
set_param_value('V_ss'  ,dynare_inputs.ss_compare(9)   );
set_param_value('d_ss'  ,dynare_inputs.ss_compare(10)   );

model; 
    % variables of interest 
    1/(C-(L^(1+(1/gamma))/(1+(1/gamma)))) = beta*(R(+1)*V(+1) + 1 -d(+1))* 1/(C(+1)-(L(+1)^(1+(1/gamma))/(1+(1/gamma))));
    L       = (W)^(gamma); 
    Y       = Z*(V*K(-1))^alpha*L^(1-alpha);
    Y       = C + I ;
    K       = (1-d)*K(-1) + I;
    Z       = Z(-1)^rho*exp(eps); 
    R       = alpha*Z*(V*K(-1))^(alpha-1)*L^(1-alpha);
    R       = delta0*omega*(V)^(omega - 1);
    W       = (1-alpha)*Z*(V*K(-1))^(alpha)*L^(-alpha);
    d       = delta0*V^(omega);
end; 

initval; 
    Z = Z_ss;
    Y = Y_ss;
    K = K_ss;
    C = C_ss;
    I = I_ss;
    L = L_ss;
    R = R_ss;
    W = W_ss;
    V = V_ss;
    d = d_ss;
end; 

shocks; 
    var eps = 1^2; 
end; 

stoch_simul(irf=20,order=1,irf_plot_threshold=0) Z R W L K I Y C V d;

