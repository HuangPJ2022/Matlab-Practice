var Z R W L K I Y C G T;
varexo T_shock; 
parameters Phi delta sigma alpha rho chi Z_ss R_ss W_ss L_ss K_ss I_ss Y_ss C_ss G_ss T_ss;

% load steady state and parameter data 
load  dynare_inputs; 

% assign parameter values 
set_param_value('Phi'  ,dynare_inputs.param_values.Phi  );
set_param_value('delta' ,dynare_inputs.param_values.delta );
set_param_value('sigma' ,dynare_inputs.param_values.sigma );
set_param_value('alpha'   ,dynare_inputs.param_values.alpha   ); 
set_param_value('rho'   ,dynare_inputs.param_values.rho   );
set_param_value('chi'   ,dynare_inputs.param_values.chi   ); 

% load steady state values as parameters
set_param_value('Z_ss'  ,dynare_inputs.ss(1)  );
set_param_value('R_ss'  ,dynare_inputs.ss(2) );
set_param_value('W_ss'  ,dynare_inputs.ss(3) );
set_param_value('L_ss'  ,dynare_inputs.ss(4)   );
set_param_value('K_ss'  ,dynare_inputs.ss(5)   ); 
set_param_value('I_ss'  ,dynare_inputs.ss(6)   );
set_param_value('Y_ss'  ,dynare_inputs.ss(7)   );
set_param_value('C_ss'  ,dynare_inputs.ss(8)   );
set_param_value('G_ss'  ,dynare_inputs.ss(9)   );
set_param_value('T_ss'  ,dynare_inputs.ss(10)   );

model; 
    1/((1+T)*C) = (1/(1+Phi)) * (1/((1+T(+1))*C(+1))) * (R(+1) + 1 - delta);
    chi/(1-L) = W / ((1+T)*C);
    Y = Z * K(-1)^alpha * L^(1-alpha);
    W = (1-alpha) * Y / L;
    R = alpha * Y / K(-1);
    K = (1-delta) * K(-1) + I;
    Y = C + I + G;
    G = T * C;
    Z = Z(-1)^(rho);
    T = T_shock; 
end;

initval; 
    Z = Z_ss;
    Y = Y_ss;
    K = K_ss;
    C = C_ss;
    I = I_ss;
    L = L_ss;
    T = T_ss;
    G = G_ss;
    R = R_ss;
    W = W_ss;
end; 


shocks; 
    var T_shock = 0.1;
end; 

stoch_simul(irf=20,order=1,irf_plot_threshold=0) Z R W L K I Y C G T;

