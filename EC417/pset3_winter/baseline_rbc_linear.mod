var c y k l i z; 
varexo eps; 
parameters beta delta eta alpha rho  I_ss Y_ss; 

% parameter calibation
beta    = 0.99; 
delta   = 0.10; 
alpha   = 0.33;
eta     = 0.25; 
rho     = 0.50;

% obtain steady state numbers from matlab 
load  dynare_inputs; 
set_param_value('Y_ss'  ,dynare_inputs.ss(2) );
set_param_value('I_ss'  ,dynare_inputs.ss(5) ); 

model(linear); 
    c(+1) -c= (1-beta*(1-delta))*(y(+1)-k); 
    c       = y-l*(1+1/eta);
    y       = z + alpha*k(-1) + (1-alpha)*l ;
    y       = c*(1-I_ss/Y_ss) + i*(I_ss/Y_ss);
    k       = (1-delta)*k(-1) + delta*i ; 
    z       = rho*z(-1) + eps; 
end; 

shocks; 
    var eps = 0.025^2; 
end; 

stoch_simul(irf=15,order=1,irf_plot_threshold=0) c y k z i l ;
close all; 


% code below is standard matlab code to produce better-looking figures 
xaxis = [0:15]; 
yaxis = ones(length(xaxis),1)*0;

subplot(2,3,1);
hold on ;
set(gca,'Box','on')
set(gca,'FontSize',16)
plot(oo_.irfs.c_eps,'-o', 'color', 'k', 'LineWidth', 1.5);
plot(xaxis,yaxis, 'color', 'r', 'LineWidth', 0.5);
hold off; 
title('Consumption', 'fontweight', 'bold');

subplot(2,3,2);
hold on ;
set(gca,'Box','on')
set(gca,'FontSize',16)
plot(oo_.irfs.y_eps,'-o', 'color', 'k', 'LineWidth', 1.5);
plot(xaxis,yaxis, 'color', 'r', 'LineWidth', 0.5);
title('Output', 'fontweight', 'bold');
hold off; 

subplot(2,3,3);
hold on ;
set(gca,'Box','on')
set(gca,'FontSize',16)
plot(oo_.irfs.k_eps,'-o', 'color', 'k', 'LineWidth', 1.5);
plot(xaxis,yaxis, 'color', 'r', 'LineWidth', 0.5);
hold off; 
title('Capital', 'fontweight', 'bold');

subplot(2,3,4);
hold on ;
set(gca,'Box','on')
set(gca,'FontSize',16)
plot(oo_.irfs.i_eps,'-o', 'color', 'k', 'LineWidth', 1.5);
plot(xaxis,yaxis, 'color', 'r', 'LineWidth', 0.5);
hold off; 
title('Investment', 'fontweight', 'bold');

subplot(2,3,5);
hold on ;
set(gca,'Box','on')
set(gca,'FontSize',16)
plot(oo_.irfs.l_eps,'-o', 'color', 'k', 'LineWidth', 1.5);
plot(xaxis,yaxis, 'color', 'r', 'LineWidth', 0.5);
hold off; 
title('Employment', 'fontweight', 'bold');
xlabel('Quarters')

subplot(2,3,6);
hold on ;
set(gca,'Box','on')
set(gca,'FontSize',16)
plot(oo_.irfs.z_eps,'-o', 'color', 'k', 'LineWidth', 1.5);
plot(xaxis,yaxis, 'color', 'r', 'LineWidth', 0.5);
hold off; 
title('Productivity', 'fontweight', 'bold');

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'baseline_RBC_linear.pdf');

