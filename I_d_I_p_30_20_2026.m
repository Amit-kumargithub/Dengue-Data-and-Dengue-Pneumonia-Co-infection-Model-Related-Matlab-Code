function I_d_I_p_30_20_2026()
%% Dengue–Pneumonia Co-infection Model with Three Control Strategies
clc; clear; close all;

%% --- Baseline Parameters
p.Lambda_h = 2500; 
p.mu_h = 1/(72*365);
p.Lambda_v = 5000; 
p.mu_v = 0.1;
p.beta_d = 0.15; 
p.beta_p = 0.08; 
p.beta_v = 0.7;
p.b = 0.5;

p.gamma_d = 1/5.5;
p.tau_d   = 1/7;
p.delta_d = 0.001;
p.omega_d = 1/180;

p.gamma_p = 1/3;
p.tau_p   = 1/14;
p.delta_p = 0.005;
p.omega_p = 1/(2*365);

p.epsilon_d = 1.5;
p.epsilon_p = 1.6;
p.epsilon_dv = 1.3;
p.rho_dp = 2.0;
p.rho_pd = 1.5;
p.phi_d = 0.5;
p.phi_p = 0.6;
p.delta_dp = 0.03;
p.tau_dp = 1/25;

p.sigma_d = 0.3; 
p.sigma_p = 0.4; 
p.sigma_v = 0.3;
p.alpha_d = 0.001; 
p.alpha_p = 0.002;

%% --- Control strategies
strategies = struct( ...
    'name', {'No Control', 'Moderate Control', 'High Control'}, ...
    'm', {0, 0.2, 0.4}, ...          % hygiene
    'c', {0, 0.15, 0.3}, ...         % contact reduction
    'kappa_v', {0, 0.2, 0.4}, ...    % vector control
    'l', {0, 0.1, 0.3});             % larval control

%% --- Initial conditions (Updated with I_d(0)=30 and I_p(0)=20)
Nh0 = 1e6;
Nv0 = 2e6;
S = Nh0 - 65;  % susceptible humans (adjusted to maintain total population)
E_d = 10; I_d = 30; R_d = 0;  % I_d(0) = 30
E_p = 10; I_p = 20; R_p = 0;  % I_p(0) = 20
I_dp = 0; V = 0;
S_v = Nv0 - 10; E_v = 5; I_v = 5;

y0 = [S E_d I_d R_d E_p I_p R_p I_dp V S_v E_v I_v];
tspan = [0 365]; % 2 years

%% --- Solve and plot
figure('Name','Control Strategy Comparison','Position',[100 100 1200 500]);

% Define colors for each strategy
colors = lines(numel(strategies)); % This creates distinct colors

for k = 1:numel(strategies)
    ctrl = strategies(k);

    % Effective transmission parameters
    p.beta_d_eff = p.beta_d * (1 - ctrl.c) * (1 - ctrl.m);
    p.beta_p_eff = p.beta_p * (1 - ctrl.c);
    p.beta_v_eff = p.beta_v * (1 - ctrl.kappa_v) * (1 - ctrl.l);
    p.Lambda_v_eff = p.Lambda_v * (1 - ctrl.l);

    % Solve ODE
    [t, y] = ode45(@(t,y) coInfectionODE(t,y,p), tspan, y0);

    % Infection curves
    I_d = y(:,3);
    I_p = y(:,6);
    I_dp = y(:,8);

    % Plot dengue with specific color
    subplot(1,2,1)
    plot(t, I_d, 'LineWidth', 2, 'Color', colors(k,:)); 
    if k == 1
        hold on;
    end

    % Plot pneumonia + co-infection with specific color
    subplot(1,2,2)
    plot(t, I_p + I_dp, 'LineWidth', 2, 'Color', colors(k,:)); 
    if k == 1
        hold on;
    end
end

% Now add titles, labels, and legends AFTER the loop
subplot(1,2,1)
xlabel('Time (days)','fontsize',14); ylabel('Dengue infected I_d','fontsize',14);
title('Dengue infection under different control strategies','fontsize',14);
legend({strategies.name},'Location','northeast','fontsize',14); 
grid on;
set(gca, 'FontSize', 14);  % Set axis tick font size

subplot(1,2,2)
xlabel('Time (days)','fontsize',14); 
ylabel('Pneumonia + Co-infection (I_p + I_{dp})','fontsize',14);
title('Pneumonia and co-infection under different control strategies','fontsize',14);
legend({strategies.name},'Location','northeast','fontsize',14); 
grid on;
set(gca, 'FontSize', 14);  % Set axis tick font size

sgtitle('Effectiveness of Control Measures (I_d(0)=30, I_p(0)=20)','fontsize',12);

end

%% --- ODE function
function dydt = coInfectionODE(~, y, p)
% Human compartments
S   = y(1);  E_d = y(2);  I_d = y(3);  R_d = y(4);
E_p = y(5);  I_p = y(6);  R_p = y(7);  I_dp = y(8);  V = y(9);
% Vector compartments
S_v = y(10); E_v = y(11); I_v = y(12);

% Total populations
N_h = S + E_d + I_d + R_d + E_p + I_p + R_p + I_dp + V;
N_v = S_v + E_v + I_v;

% Forces of infection
lambda_d = p.beta_d_eff * (I_d + p.epsilon_d*I_dp)/N_h;
lambda_p = p.beta_p_eff * (I_p + p.epsilon_p*I_dp)/N_h;
lambda_v = p.beta_v_eff * p.b * (I_d + p.epsilon_dv*I_dp)/N_h;

% Human dynamics
dS = p.Lambda_h - (p.beta_d*I_v*S)/N_h - lambda_p*S - (p.mu_h + p.alpha_d + p.alpha_p)*S;

dE_d = (p.beta_d*I_v*(S + p.sigma_p*R_p + p.sigma_v*V))/N_h - (p.gamma_d + p.mu_h)*E_d ...
       - p.rho_dp*lambda_p*E_d;

dI_d = p.gamma_d*E_d - (p.mu_h + p.delta_d + p.tau_d)*I_d ...
       - p.rho_dp*lambda_p*I_d;

dR_d = p.tau_d*I_d - (p.mu_h + p.omega_d)*R_d ...
       - p.phi_d*lambda_p*R_d;

dE_p = lambda_p*(S + p.sigma_d*R_d + p.sigma_v*V) - (p.gamma_p + p.mu_h)*E_p ...
       - (p.rho_pd*p.beta_d*I_v*E_p)/N_h;

dI_p = p.gamma_p*E_p - (p.mu_h + p.delta_p + p.tau_p)*I_p ...
       - (p.rho_pd*p.beta_d*I_v*I_p)/N_h;

dR_p = p.tau_p*I_p - (p.mu_h + p.omega_p)*R_p ...
       - (p.phi_p*p.beta_d*I_v*R_p)/N_h;

dI_dp = p.rho_dp*lambda_p*(E_d + I_d) + (p.rho_pd*p.beta_d*I_v*(E_p + I_p))/N_h ...
      + p.phi_d*lambda_p*R_d +( p.phi_p*p.beta_d*I_v*R_p)/N_h ...
      - (p.mu_h + p.delta_dp + p.tau_dp)*I_dp;

dV = p.alpha_d*S + p.alpha_p*S - p.mu_h*V ...
     - p.sigma_v* lambda_p*V-(p.sigma_v*p.beta_d*I_v*V)/N_h;

% Mosquito dynamics
dS_v = p.Lambda_v_eff - lambda_v*S_v - p.mu_v*S_v;
dE_v = lambda_v*S_v - (p.gamma_d + p.mu_v)*E_v;
dI_v = p.gamma_d*E_v - p.mu_v*I_v;

dydt = [dS; dE_d; dI_d; dR_d; dE_p; dI_p; dR_p; dI_dp; dV; dS_v; dE_v; dI_v];
end