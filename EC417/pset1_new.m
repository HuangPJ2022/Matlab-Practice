% Author: Louie Po-Jui Huang
% Date: 13 Oct 2025
% Des: This script is for EC417 pset1.

%% Parameters
beta  = 0.96;     % Subjective discount factor
delta = 0.08;     % Depreciation rate

%% Compute the optimal and maximal steady state capital(Qa, Qb)
k_ss = (3*((1/beta)-1+delta))^(-3/2); % optimal ss capital(Qa)
k_G = (1/delta)^(3/2); % maximal ss capital(Qb)

%% Generate initial consumption vector
k_0 = 0.1*k_ss; % given initial capital
c_max = (k_0)^(1/3) + (1-delta)*k_0; % compute maximal consumption
temp_vec = linspace(0, c_max, 22); % pretend this is an open interval
c0_vec = temp_vec(2 : end - 1); % all values but both sides

%% Simulation Model
% Prepare container
k_paths = zeros(100, 20);
c_paths = zeros(100, 20);

% Iterate
for i = 1:20
    c_paths(1, i) = c0_vec(i);
    k_paths(1, i) = k_0;
    
    for t = 1:99
        % Stop simulating this path if capital <0
        if k_paths(t, i) < 0
            fprintf('%d-th path failed', i)
            fprintf('at %d period.\n', t)
            break;
        end

        % Capital constraint
        k_paths(t+1, i) = k_paths(t, i)^(1/3) + (1-delta)*k_paths(t, i) - c_paths(t, i);
        % EE
        c_paths(t+1, i) = c_paths(t, i) * beta * ((1/3)*(k_paths(t+1, i))^(-2/3) + 1 - delta);
    end
end

% Plot
% Prepare legend label
for i = 1:7
    legend_labels{i} = ['Path = ', num2str(i)];
end

% Create a figure and set subplots
figure;

subplot(1, 2, 1); % subplots1
plot(1:100, c_paths(:,1:7)); % plot consumption paths
hold on;
title('Consumption Paths');
xlabel('Time');
ylabel('Consumption (c)');
legend(legend_labels, 'Location', 'best');
grid on;

subplot(1, 2, 2); % subplots2
plot(1:100, k_paths(:,1:7)); % plot capital paths
hold on;
title('Capital Paths');
xlabel('Time');
ylabel('Capital (k)');
legend(legend_labels, 'Location', 'best');
grid on;

save pset1_4d.png

%% Compute TVC values
% Prepare container
tvc_paths = zeros(99,8);

% Iterate
for i = 1:7
    for t = 1: 99
        tvc_paths(t,i) = (beta)^(t) * k_paths(t + 1,i) * (c_paths(t, i))^(-1);
    end
end

% Plot
figure;
plot(1:99, tvc_paths(:,1:7)); % plot TVC paths
title('TVC Paths');
xlabel('Time');
ylabel('TVC');
legend(legend_labels, 'Location', 'best');
grid on;

save pset1_4e.png

