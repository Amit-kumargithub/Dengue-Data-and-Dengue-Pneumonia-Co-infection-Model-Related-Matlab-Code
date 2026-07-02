%% DENGUE–PNEUMONIA CO-INFECTION MODEL WITH SEASONALITY
% Periodic epidemic scenario driven by mosquito biting seasonality
% Simulates 6 years (2190 days)
% Author: ChatGPT (GPT-5)
function seasonal_coinfection_model

clear; clc; close all;

%% ===================== PARAMETERS =========================
p.Lambda_h = 1000;             % human recruitment rate
p.mu_h = 1/(72*365);           % natural mortality (humans)
p.Lambda_v = 15000;            % mosquito recruitment
p.mu_v = 0.08;                 % mosquito mortality
p.beta_d = 0.25;               % dengue transmission (moderate)
p.beta_p = 0.07;               % pneumonia transmission
p.beta_v = 0.8;                % vector-to-human transmission
p.b0 = 0.6;                    % average mosquito biting rate
p.eta = 0.35;                  % seasonality amplitude (35%)
p.gamma_d = 1/6;               % incubation (dengue)
p.tau_d = 1/10;                % recovery (dengue)
p.delta_d = 0.0008;            % dengue death
p.omega_d = 1/180;             % waning (dengue)
p.gamma_p = 1/3;               % incubation (pneumonia)
p.tau_p = 1/14;                % recovery (pneumonia)
p.delta_p = 0.004;             % pneumonia death
p.omega_p = 1/(1.5*365);       % waning (pneumonia)
p.epsilon_d = 1.5;             % co-infection enhancement dengue
p.epsilon_p = 1.4;             % co-infection enhancement pneumonia
p.epsilon_dv = 1.3;            % mosquito enhancement
p.rho_dp = 1.8;                % dengue→pneumonia susceptibility
p.rho_pd = 1.5;                % pneumonia→dengue susceptibility
p.phi_d = 0.6; p.phi_p = 0.5;  % reinfection risk
p.delta_dp = 0.02;             % co-infection death
p.tau_dp = 1/25;               % recovery co-infection
p.sigma_d = 0.3; p.sigma_p = 0.4; p.sigma_v = 0.3;
p.alpha_d = 0.0005; p.alpha_p = 0.001;  % mild vaccination

%% ===================== INITIAL CONDITIONS =========================
Nh0 = 1e6; Nv0 = 1.2e6;    % initial populations
S0   = Nh0 - (200+100+50+20+0+0+0);
Ed0  = 100; Id0 = 200; Rd0 = 0;
Ep0  = 50;  Ip0 = 100; Rp0 = 0;
Idp0 = 20;  V0  = 0;
Sv0 = Nv0 - (400+200);
Ev0 = 400; Iv0 = 200;
y0 = [S0 Ed0 Id0 Rd0 Ep0 Ip0 Rp0 Idp0 V0 Sv0 Ev0 Iv0];

%% ===================== SIMULATION =========================
tspan = [0 3*365];  % simulate 6 years
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

[t, y] = ode45(@(t,y) rhs_periodic(t,y,p), tspan, y0, opts);

% Extract solutions
S=y(:,1); Ed=y(:,2); Id=y(:,3); Rd=y(:,4);
Ep=y(:,5); Ip=y(:,6); Rp=y(:,7); Idp=y(:,8); V=y(:,9);
Sv=y(:,10); Ev=y(:,11); Iv=y(:,12);
Nh = S+Ed+Id+Rd+Ep+Ip+Rp+Idp+V;

%% ===================== PLOTS =========================
figure('Position',[100 100 1200 750])

subplot(2,2,1)
plot(t/365, Id./Nh, 'r', 'LineWidth',1.8); hold on
plot(t/365, Ip./Nh, 'b', 'LineWidth',1.8);
plot(t/365, Idp./Nh, 'm', 'LineWidth',1.5);
xlabel('Time (years)','FontSize',14)
ylabel('Prevalence (fraction of humans)','FontSize',14)
legend('Dengue (I_d)','Pneumonia (I_p)','Co-infected (I_{dp})','Location','best','FontSize',14)
title('Oscillatory Infection Dynamics','FontWeight','bold','FontSize',14)
set(gca, 'FontSize', 13, 'LineWidth', 2);

subplot(2,2,2)
plot(t/365, Iv./(Sv+Ev+Iv), 'k','LineWidth',1.8)
xlabel('Time (years)','FontSize',14)
ylabel('Infected Mosquito Fraction','FontSize',14)
title('Vector Infection Cycles','FontWeight','bold','FontSize',14)
set(gca, 'FontSize', 13, 'LineWidth', 2);

subplot(2,2,3)
plot(t/365, S./Nh, 'g', t/365, (Rd+Rp)./Nh, 'c', 'LineWidth',1.5)
xlabel('Time (years)','FontSize',14)
ylabel('Fraction of Humans','FontSize',14)
legend('Susceptible','Recovered','Location','best','FontSize',14)
title('Host Immunity and Susceptibility','FontWeight','bold','FontSize',14)
set(gca, 'FontSize', 13, 'LineWidth', 2);

subplot(2,2,4)
b_t = p.b0*(1 + p.eta*sin(2*pi*t/365));
plot(t/365, b_t, 'k','LineWidth',1.8)
xlabel('Time (years)','FontSize',14)
ylabel('Biting rate b(t)','FontSize',14)
title('Seasonal Forcing','FontWeight','bold','FontSize',14)
set(gca, 'FontSize', 13, 'LineWidth', 2);

%title('Seasonally Forced Dengue–Pneumonia Co-infection Model (Periodic Scenario)','FontWeight','bold')

%% ===================== NESTED FUNCTION =========================
function dydt = rhs_periodic(t, y, p)
% Differential equations with seasonal biting rate

S=y(1); Ed=y(2); Id=y(3); Rd=y(4);
Ep=y(5); Ip=y(6); Rp=y(7); Idp=y(8); V=y(9);
Sv=y(10); Ev=y(11); Iv=y(12);

Nh = S+Ed+Id+Rd+Ep+Ip+Rp+Idp+V;

% Seasonal biting rate
b = p.b0*(1 + p.eta*sin(2*pi*t/365));

% Forces of infection
lambda_d = p.beta_d * (Id + p.epsilon_d*Idp)/Nh;
lambda_p = p.beta_p * (Ip + p.epsilon_p*Idp)/Nh;
lambda_v = p.beta_v * b * (Id + p.epsilon_dv*Idp)/Nh;

% Human dynamics
dS = p.Lambda_h - (p.beta_d*Iv*S)/Nh - lambda_p*S ...
     - (p.mu_h + p.alpha_d + p.alpha_p)*S; ...;

dEd = (p.beta_d*Iv*(S + p.sigma_p*Rp + p.sigma_v*V))/Nh ...
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
     - p.rho_pd*Ep*(p.beta_d*Iv*S)/Nh;

dIp = p.gamma_p*Ep ...
     - (p.mu_h + p.delta_p + p.tau_p)*Ip ...
     - p.rho_pd*Ip*(p.beta_d*Iv*S)/Nh;

dRp = p.tau_p*Ip ...
     - (p.mu_h)*Rp ...
     - p.phi_p*Rp*(p.beta_d*Iv*S)/Nh;

dIdp = p.rho_dp*lambda_p*(Ed + Id) ...
      + p.rho_pd*(Ep + Ip)*(p.beta_d*Iv*S)/Nh ...
      + p.phi_d*lambda_p*Rd + p.phi_p*Rp*(p.beta_d*Iv*S)/Nh ...
      - (p.mu_h + p.delta_dp + p.tau_dp)*Idp;

dV = (p.alpha_d + p.alpha_p)*S ...
    - p.mu_h*V ...
    - p.sigma_v*((p.beta_d*Iv*S)/Nh + lambda_p)*V;

% Mosquito dynamics
dSv = p.Lambda_v - lambda_v*Sv - p.mu_v*Sv;
dEv = lambda_v*Sv - (p.gamma_d + p.mu_v)*Ev;
dIv = p.gamma_d*Ev - p.mu_v*Iv;

dydt = [dS; dEd; dId; dRd; dEp; dIp; dRp; dIdp; dV; dSv; dEv; dIv];