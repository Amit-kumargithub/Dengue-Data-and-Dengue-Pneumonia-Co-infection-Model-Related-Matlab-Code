%% DENGUE–PNEUMONIA CO-INFECTION MODEL (OUTBREAK SCENARIO)
% Self-contained MATLAB script
% Simulates 1 years (365 days)
% a single pick version
function dengue_pneumonia_coinfection_time_dynamics
clear; clc; close all;

%% ===================== PARAMETERS =========================
p.Lambda_h = 2000;              % human recruitment rate
p.mu_h = 1/(72*365);            % natural mortality (humans)
p.Lambda_v = 20000;             % mosquito recruitment rate (↑)
p.mu_v = 0.1;                   % mosquito mortality
p.beta_d = 0.35;                % dengue transmission (↑)
p.beta_p = 0.10;                % pneumonia transmission
p.beta_v = 0.9;                 % vector-to-human transmission (↑)
p.b = 0.8;                      % mosquito biting rate (↑)
p.gamma_d = 1/5.5;              % incubation rate (dengue)
p.tau_d = 1/14;                 % recovery (slower)
p.delta_d = 0.001;              % dengue-induced death
p.omega_d = 1/180;              % waning (dengue)
p.gamma_p = 1/3;                % incubation (pneumonia)
p.tau_p = 1/21;                 % recovery (slower)
p.delta_p = 0.005;              % pneumonia death
p.omega_p = 1/(2*365);          % waning (pneumonia)
p.epsilon_d = 1.6;              % dengue infectiousness enhancement
p.epsilon_p = 1.6;              % pneumonia infectiousness enhancement
p.epsilon_dv = 1.4;             % mosquito infection enhancement
p.rho_dp = 2.0;                 % dengue→pneumonia susceptibility
p.rho_pd = 1.8;                 % pneumonia→dengue susceptibility
p.phi_d = 0.6;                  % reinfection risk (dengue)
p.phi_p = 0.6;                  % reinfection risk (pneumonia)
p.delta_dp = 0.03;              % co-infection death
p.tau_dp = 1/25;                % recovery from co-infection
p.sigma_d = 0.3; p.sigma_p = 0.4; p.sigma_v = 0.3;
p.alpha_d = 0.001; p.alpha_p = 0.002;  % vaccination (low)
%p.alpha_d = 0.005; p.alpha_p = 0.008;  % vaccination (high)
p.kappa_v = 0.3; p.l = 0.2; p.m = 0.2; p.c = 0.15;

%% ===================== INITIAL CONDITIONS =========================
Nh0 = 1e6;     % human population
Nv0 = 2e6;     % mosquito population (↑)

S0   = Nh0 - (500+300+200+50+0+0+0); 
Ed0  = 200;  Id0 = 300;  Rd0 = 0;
Ep0  = 100;  Ip0 = 200;  Rp0 = 0;
Idp0 = 50;   V0  = 0;

Sv0 = Nv0 - (1000+500);
Ev0 = 1000; Iv0 = 500;

y0 = [S0 Ed0 Id0 Rd0 Ep0 Ip0 Rp0 Idp0 V0 Sv0 Ev0 Iv0];

%% ===================== SIMULATION =========================
tspan = [0 1*400];  % days
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

[t, y] = ode45(@(t,y) dengue_pneumonia_rhs(t,y,p), tspan, y0, opts);

S=y(:,1); Ed=y(:,2); Id=y(:,3); Rd=y(:,4);
Ep=y(:,5); Ip=y(:,6); Rp=y(:,7); Idp=y(:,8); V=y(:,9);
Sv=y(:,10); Ev=y(:,11); Iv=y(:,12);
Nh = S+Ed+Id+Rd+Ep+Ip+Rp+Idp+V;

%% ===================== PLOTS =========================
figure('Position',[100 100 1200 700])

subplot(2,2,1)
plot(t, [S Id Ip Idp]/1e5, 'LineWidth',1.8)
xlabel('Days','FontSize',11)
ylabel('Humans (\times10^5)','FontSize',14)
legend('S','I_d','I_p','I_{dp}','Location','best')
title('Human Infection Dynamics','FontWeight','bold','FontSize',14)
set(gca, 'FontSize', 13, 'LineWidth', 2);

subplot(2,2,2)
plot(t, [Ed Ep Idp], 'LineWidth',1.8)
xlabel('Days','FontSize',11)
ylabel('Individuals','FontSize',14)
legend('E_d','E_p','I_{dp}','Location','best')
title('Exposed and Co-infected Dynamics','FontWeight','bold','FontSize',14)
set(gca, 'FontSize', 13, 'LineWidth', 2);

subplot(2,2,3)
plot(t, [Rd Rp V], 'LineWidth',1.8)
xlabel('Days','FontSize',11)
ylabel('Individuals','FontSize',14)
legend('R_d','R_p','V','Location','best')
title('Recovered and Vaccinated','FontWeight','bold','FontSize',14)
set(gca, 'FontSize', 13, 'LineWidth', 2);

subplot(2,2,4)
plot(t, [Sv Iv]/1e5, 'LineWidth',1.8)
xlabel('Days','FontSize',11)
ylabel('Mosquito Pop. (\times10^5)','FontSize',14)
legend('S_v','I_v','Location','best')

title('Mosquito Population Dynamics','FontWeight','bold','FontSize',14)
set(gca, 'FontSize', 13, 'LineWidth', 2);
%title('Large Outbreak Scenario: Dengue–Pneumonia Co-infection (3 Years)','FontWeight','bold')

%% ===================== NESTED FUNCTION =========================
end
function dydt = dengue_pneumonia_rhs(~, y, p)
% Right-hand side of the dengue–pneumonia co-infection system

S=y(1); Ed=y(2); Id=y(3); Rd=y(4);
Ep=y(5); Ip=y(6); Rp=y(7); Idp=y(8); V=y(9);
Sv=y(10); Ev=y(11); Iv=y(12);

Nh = S+Ed+Id+Rd+Ep+Ip+Rp+Idp+V;

% Forces of infection (humans)
lambda_d = p.beta_d * Iv/Nh;
lambda_p = p.beta_p * (Ip + p.epsilon_p*Idp)/Nh;

% Humans
dS = p.Lambda_h - lambda_d*S - lambda_p*S ...
     - (p.mu_h + p.alpha_d + p.alpha_p)*S ...;

dEd = lambda_d*(S + p.sigma_p*Rp + p.sigma_v*V) ...
     - (p.gamma_d + p.mu_h)*Ed ...
     - p.rho_dp*lambda_p*Ed;

dId = p.gamma_d*Ed ...
     - (p.mu_h + p.delta_d + p.tau_d)*Id ...
     - p.rho_dp*lambda_p*Id;

dRd = p.tau_d*Id ...
     - (p.mu_h)*Rd ...
     - p.phi_d*lambda_p*Rd;

dEp = lambda_p*(S + p.sigma_d*Rd + p.sigma_v*V) ...
     - (p.gamma_p + p.mu_h)*Ep ...
     - p.rho_pd*lambda_d*Ep;

dIp = p.gamma_p*Ep ...
     - (p.mu_h + p.delta_p + p.tau_p)*Ip ...
     - p.rho_pd*lambda_d*Ip;

dRp = p.tau_p*Ip ...
     - (p.mu_h)*Rp ...
     - p.phi_p*lambda_d*Rp;

dIdp = p.rho_dp*lambda_p*(Ed + Id) ...
      + p.rho_pd*lambda_d*(Ep + Ip) ...
      + p.phi_d*lambda_p*Rd + p.phi_p*lambda_d*Rp ...
      - (p.mu_h + p.delta_dp + p.tau_dp)*Idp;

dV = (p.alpha_d + p.alpha_p)*S ...
    - p.mu_h*V ...
    - p.sigma_v*(lambda_d + lambda_p)*V;

% Mosquito dynamics
lambda_v = p.beta_v*p.b*(Id + p.epsilon_dv*Idp)/Nh;
dSv = p.Lambda_v - lambda_v*Sv - p.mu_v*Sv;
dEv = lambda_v*Sv - (p.gamma_d + p.mu_v)*Ev;
dIv = p.gamma_d*Ev - p.mu_v*Iv;

dydt = [dS; dEd; dId; dRd; dEp; dIp; dRp; dIdp; dV; dSv; dEv; dIv];
end