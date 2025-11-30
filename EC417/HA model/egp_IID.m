% Endogenous Grid Points with IID Income
% Greg Kaplan 2017

clear;
close all;

%% PARAMETERS

% preferences
risk_aver   = 2;
beta        = 0.961;

%returns
r           = 0.04;
R = 1+ r;

% income risk: discretized N(mu,sigma^2)
mu_y    = 1;
sd_y    = 0.2;
ny      = 5;

% asset gridss
na          = 500;
amax        = 50; 
borrow_lim  = 0;
agrid_par   = 0.5; %1 for linear, 0 for L-shaped

% computation
max_iter    = 1000;
tol_iter    = 1.0e-6;
Nsim        = 50000;
Tsim        = 1000;

%mpc options
mpcamount1  = 1.0e-10; %approximate thoeretical mpc
mpcamount2  = 0.10; % one percent of average income: approx $500

%% OPTIONS
Display     = 1;
DoSimulate  = 1;
MakePlots   = 1;
ComputeMPC  = 1;


%% DRAW RANDOM NUMBERS
rng(2017);
yrand = rand(Nsim,Tsim);

%% SET UP GRIDS

% assets
agrid = linspace(0,1,na)';
agrid = agrid.^(1./agrid_par);
agrid = borrow_lim + (amax-borrow_lim).*agrid;

% income: disretize normal distribution
width = fzero(@(x)discrete_normal(ny,mu_y,sd_y,x),2);
[temp,ygrid,ydist] = discrete_normal(ny,mu_y,sd_y,width);
ycumdist = cumsum(ydist);

%% UTILITY FUNCTION

u = @(c)(c.^(1-risk_aver)-1)./(1-risk_aver);
u1 = @(c) c.^(-risk_aver); %mu
u1inv = @(u) u.^(-1./risk_aver); %inverse function, plug in mu to get c

%% INITIALIZE CONSUMPTION FUNCTION

conguess = zeros(na,ny);
for iy = 1:ny
    conguess(:,iy) = r.*agrid+ygrid(iy);
end


%% ITERATE ON EULER EQUATION WITH ENDOGENOUS GRID POINTS

con = conguess;

iter = 0;
cdiff = 1000;

while iter <= max_iter && cdiff>tol_iter
    iter = iter + 1;
    sav = zeros(na,ny);
    
    conlast = con;
    
    emuc = u1(conlast)*ydist; %expected mu
    muc1 = beta.*R.*emuc; %EE
    con1 = u1inv(muc1); % plug in inverse function to get consumption
    
    % loop over income
    for iy = 1:ny
        
        ass1(:,iy) = (con1 + agrid -ygrid(iy))./R; %BC to get asset tomorrow
        
        % loop over current period ssets
        for ia  = 1:na 
            if agrid(ia)<ass1(1,iy) %borrowing constraint binds
                sav(ia,iy) = borrow_lim;
                
            else %borrowing constraint does not bind;
                sav(ia,iy) = lininterp1(ass1(:,iy),agrid,agrid(ia));
                
            end                
        end
        con(:,iy) = R.*agrid +ygrid(iy) - sav(:,iy);
    end   
    
    cdiff = max(max(abs(con-conlast)));
    disp(['Iteration no. ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
end    


%% SIMULATE
if DoSimulate==1

    yindsim = zeros(Nsim,Tsim);
    asim = zeros(Nsim,Tsim);
    
    %create interpolating function
    for iy = 1:ny
        savinterp{iy} = griddedInterpolant(agrid,sav(:,iy),'linear');
    end
    
    %loop over time periods
    for it = 1:Tsim
        if Display >=1 && mod(it,100) ==0
            disp([' Simulating, time period ' int2str(it)]);
        end
        
        %income realization: note we vectorize simulations at once because
        %of matlab, in other languages we would loop over individuals
        yindsim(yrand(:,it)<= ycumdist(1),it) = 1;
        for iy = 2:ny
            yindsim(yrand(:,it)> ycumdist(iy-1) & yrand(:,it)<=ycumdist(iy),it) = iy;
        end
        
        % asset choice
        if it<Tsim
            for iy = 1:ny
                asim(yindsim(:,it)==iy,it+1) = savinterp{iy}(asim(yindsim(:,it)==iy,it));
            end
        end
    end
    
    %assign actual income values;
    ysim = ygrid(yindsim);

end

 
%% MAKE PLOTS
if MakePlots ==1 
    figure(1);
    
    % consumption policy function
    subplot(2,4,1);
    plot(agrid,con(:,1),'b-',agrid,con(:,ny),'r-','LineWidth',1);
    grid;
    xlim([0 amax]);
    title('Consumption Policy Function');
    legend('Lowest income state','Highest income state');

    % savings policy function
    subplot(2,4,2);
    plot(agrid,sav(:,1)-agrid,'b-',agrid,sav(:,ny)-agrid,'r-','LineWidth',1);
    hold on;
    plot(agrid,zeros(na,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim([0 amax]);
    title('Savings Policy Function (a''-a)');
    
    % consumption policy function: zoomed in
    subplot(2,4,3);
    plot(agrid,con(:,1),'b-o',agrid,con(:,ny),'r-o','LineWidth',2);
    grid;
    xlim([0 1]);
    title('Consumption: Zoomed');
    
     % savings policy function: zoomed in
    subplot(2,4,4);
    plot(agrid,sav(:,1)-agrid,'b-o',agrid,sav(:,ny)-agrid,'r-o','LineWidth',2);
    hold on;
    plot(agrid,zeros(na,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim([0 1]);
    title('Savings: Zoomed (a''-a)');
    
    
    %income distribution
    subplot(2,4,5);
    hist(ysim(:,Tsim),ygrid);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[0 0.5 0.5],'EdgeColor','blue','LineStyle','-');
    ylabel('')
    title('Income distribution');
    
    %asset distribution
    subplot(2,4,6:7);
    hist(asim(:,Tsim),100);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[.7 .7 .7],'EdgeColor','black','LineStyle','-');
    ylabel('')
    title('Asset distribution');

    %convergence check
    subplot(2,4,8);
    plot([1:Tsim]',mean(asim,1),'k-','LineWidth',1.5);
    xlabel('Time Period');
    title('Mean Asset Convergence');
    
    
   % asset distribution statistics
    aysim = asim(:,Tsim) ./ mean(ysim(:,Tsim));
    disp(['Mean assets  (relative to mean income) : ' num2str(mean(aysim))]);
    disp(['Fraction borrowing constrained: ' num2str(sum(aysim==borrow_lim)./Nsim * 100) '%']);
    disp(['10th Percentile: ' num2str(quantile(aysim,.1))]);
    disp(['50th Percentile: ' num2str(quantile(aysim,.5))]);
    disp(['90th Percentile: ' num2str(quantile(aysim,.9))]);
    disp(['99th Percentile: ' num2str(quantile(aysim,.99))]);

    %table of parameters;
    tabparam  = zeros(2,1);
    tabparam(1) = beta; %discount factor
    tabparam(2) = risk_aver; %risk aversion

    %table of statistics
    tabstat = zeros(12,1);
    tabstat(1) = var(log(ysim(:,Tsim))); %variance log earnings;
    tabstat(2) = ginicoeff(ysim(:,Tsim)); %gini earnings;
    tabstat(3) = mean(aysim); %mean wealth (relative to mean earnings);
    tabstat(4) = quantile(aysim,.5); %median wealth (relative to mean earnings);
    tabstat(5) = ginicoeff(aysim); %wealth gini;
    tabstat(6) = quantile(aysim,.9) ./ quantile(aysim,.5); %90-10 wealth distribution
    tabstat(7) = quantile(aysim,.99) ./ quantile(aysim,.5); %99-50 wealth distribution
    tabstat(8) = sum(aysim<=0)./Nsim; %fraction  less than equal to zero;
    tabstat(9) = sum(aysim<=0.05)./Nsim; %fraction  less than equal to 5% av. annual gross earnings;
    tabstat(10) = sum(aysim(aysim>=quantile(aysim,.9))) ./ sum(aysim); %top 10% wealth share;
    tabstat(11) = sum(aysim(aysim>=quantile(aysim,.99))) ./ sum(aysim); %top 1% wealth share;
    tabstat(12) = sum(aysim(aysim>=quantile(aysim,.999))) ./ sum(aysim); %top 0.1% wealth share;
    
    taboutput = [tabparam; tabstat];

end

%% COMPUTE MPCs
if ComputeMPC ==1
    
    %theoretical mpc lower bound
    mpclim = R*((beta*R)^-(1./risk_aver))-1;
    
    for iy = 1:ny
        %create interpolating function
        coninterp{iy} = griddedInterpolant(agrid,con(:,iy),'linear');
        
        mpc1(:,iy) = ( coninterp{iy}(agrid+mpcamount1) - con(:,iy) ) ./ mpcamount1;
        mpc2(:,iy) = ( coninterp{iy}(agrid+mpcamount2) - con(:,iy) ) ./ mpcamount2;
        
    end

    
    figure(2);
    
    % mpc functions
    subplot(1,2,1);
    plot(agrid,mpc1(:,1),'b-',agrid,mpc2(:,1),'b--',...
        agrid,mpc1(:,ny),'r-',agrid,mpc2(:,ny),'r--','LineWidth',1);
    hold on;
    plot(agrid,mpclim.*ones(size(agrid)),'k:','LineWidth',2);    
    hold off;
    grid on;
    xlim([0 10]);
    title('MPC Function');
    legend('Lowest income state: amount 1','Lowest income state: amount 2',...
        'Highest income state: amount 1','Highest income state: amount 2',...
        ['Theoretical MPC limit = ', num2str(mpclim)],'Location','NorthEast');
    
    % mpc distribution
    mpc1sim = zeros(Nsim,1);
    mpc2sim = zeros(Nsim,1);
    for iy = 1:ny
        mpc1sim(yindsim(:,Tsim)==iy) = ( coninterp{iy}(asim(yindsim(:,Tsim)==iy,Tsim)+mpcamount1) - coninterp{iy}(asim(yindsim(:,Tsim)==iy,Tsim)) ) ./ mpcamount1;
        mpc2sim(yindsim(:,Tsim)==iy) = ( coninterp{iy}(asim(yindsim(:,Tsim)==iy,Tsim)+mpcamount2) - coninterp{iy}(asim(yindsim(:,Tsim)==iy,Tsim)) ) ./ mpcamount2;
    end
    
    subplot(1,2,2);
    histogram(mpc1sim,[0:0.02:1.5]);
    grid on;
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[.7 .7 .7],'EdgeColor','black','LineStyle','-');
    ylabel('')
    title('MPC distribution');
    
    % mpc distribution statistics
    disp(['Mean MPC amount 1: ' num2str(mean(mpc1sim))]);
    disp(['Mean MPC amount 2: ' num2str(mean(mpc2sim))]);

    
end    
    
    
%% Lorenz curve
% 1. Sort the asset data
A_sorted = sort(aysim);

% 2. Calculate Cumulative Population and Wealth
pop_share = (1:length(A_sorted))' ./ length(A_sorted);
wealth_share = cumsum(A_sorted) ./ sum(A_sorted);

% 3. Plot
figure;
plot(pop_share, wealth_share, 'LineWidth', 2); hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1); % Line of Equality
xlabel('Cumulative Share of Population');
ylabel('Cumulative Share of Wealth');
title('Lorenz Curve of Wealth');
legend('Model', 'Perfect Equality', 'Location', 'NorthWest');
grid on;

%% (d)
% --- 1. Filter and Sort the Model Data ---
% We only want the Top 10%, as requested
cutoff_90 = quantile(aysim, 0.90);
a_tail = aysim(aysim >= cutoff_90);
a_tail = sort(a_tail); % Sort ascending

% --- 2. Calculate the Survival Probability (1 - G(a)) ---
% The probability of having wealth > a_i
N_tail = length(a_tail);
% Create a ranking from 1 (lowest in tail) to N (highest)
rank = (1:N_tail)';
% Probability of being ABOVE this rank
% If you are the richest person, probability is 1/N. 
% If you are the poorest in the tail, probability is near 1 (relative to the tail)
survival_prob = 1 - (rank / (N_tail + 1)); 

% Normalize to the total population share (since this is the top 10%)
% The survival prob relative to the WHOLE population ranges from 0.1 to 0
survival_prob_total = survival_prob * 0.10;

% --- 3. Log-Log Transformation ---
log_a_model = log(a_tail);
log_survival_model = log(survival_prob_total);

% --- 4. Plotting ---
figure;
hold on;

% Plot YOUR MODEL (Blue)
plot(log_a_model, log_survival_model, 'b-', 'LineWidth', 2);

% Plot THE DATA from Problem Set 6 (Item d)
% You likely have vectors like 'wealth_data_sorted' or similar from the previous prob set.
% Assuming you have variables: log_a_data, log_survival_data
% plot(log_a_data, log_survival_data, 'r--', 'LineWidth', 2); 

xlabel('log(Wealth)');
ylabel('log(1 - G(a))');
title('Tail of Wealth Distribution: Model vs. Data');
legend('Aiyagari Model', 'US Data (Power Law)');
grid on;
hold off;

    
    
    
    
    
    
    
