%% Dengue-Pneumonia Optimal Control - WITH THREE SCENARIOS
clear; clc; close all;

%% Parameters
Lambda_h = 2000; mu_h = 1/(72*365);
gamma_d = 1/5.5; tau_d = 1/7; delta_d = 0.001;
gamma_p = 1/3; tau_p = 1/14; delta_p = 0.005;
rho_dp = 2; rho_pd = 1.5; phi_d = 0.5; phi_p = 0.6;
epsilon_p = 1.6; epsilon_dv = 1.3; sigma_v = 0.3;
alpha_d = 0.001; alpha_p = 0.002; delta_dp = 0.01; tau_dp = 1/10;
Lambda_v = 5000; mu_v = 0.1; gamma_v = 1/10; b = 0.5;
beta_d = 0.15; beta_p = 0.08; beta_v = 0.7;
sigma_d = 0.3; sigma_p = 0.4;

% Weights for optimal control
A1 = 10; A2 = 10; A3 = 10;
B1 = 1; B2 = 1; B3 = 1; B4 = 1;

% Initial conditions
y0 = [950; 0; 300; 0; 0; 20; 0; 0; 0; 4970; 0; 30];

%% Time
T = 250; N = 500; t = linspace(0,T,N); dt = t(2)-t(1);

%% ==================== SCENARIO 1: u = 0 (No control) ====================
fprintf('Simulating scenario: u = 0...\n');
Y_u0 = zeros(N,12);
Y_u0(1,:) = y0';

for i = 1:N-1
    y = Y_u0(i,:)';
    S=y(1); Ed=y(2); Id=y(3); Rd=y(4); Ep=y(5); Ip=y(6); Rp=y(7); Idp=y(8); V=y(9); Sv=y(10); Ev=y(11); Iv=y(12);
    Nh = max(sum(y(1:9)),1);
    
    FOI_d = beta_d * Iv / Nh;
    FOI_p = beta_p * (Ip + epsilon_p*Idp) / Nh;
    FOI_v = beta_v * b * (Id + epsilon_dv*Idp) / Nh;
    
    k1 = zeros(12,1);
    k1(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
    k1(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
    k1(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
    k1(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
    k1(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
    k1(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
    k1(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
    k1(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
    k1(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
    k1(10) = Lambda_v - FOI_v*Sv - mu_v*Sv;
    k1(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
    k1(12) = gamma_v*Ev - mu_v*Iv;
    
    y2 = y + dt/2*k1;
    S=y2(1); Ed=y2(2); Id=y2(3); Rd=y2(4); Ep=y2(5); Ip=y2(6); Rp=y2(7); Idp=y2(8); V=y2(9); Sv=y2(10); Ev=y2(11); Iv=y2(12);
    Nh = max(sum(y2(1:9)),1);
    FOI_d = beta_d * Iv / Nh;
    FOI_p = beta_p * (Ip + epsilon_p*Idp) / Nh;
    FOI_v = beta_v * b * (Id + epsilon_dv*Idp) / Nh;
    
    k2 = zeros(12,1);
    k2(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
    k2(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
    k2(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
    k2(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
    k2(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
    k2(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
    k2(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
    k2(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
    k2(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
    k2(10) = Lambda_v - FOI_v*Sv - mu_v*Sv;
    k2(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
    k2(12) = gamma_v*Ev - mu_v*Iv;
    
    y3 = y + dt/2*k2;
    S=y3(1); Ed=y3(2); Id=y3(3); Rd=y3(4); Ep=y3(5); Ip=y3(6); Rp=y3(7); Idp=y3(8); V=y3(9); Sv=y3(10); Ev=y3(11); Iv=y3(12);
    Nh = max(sum(y3(1:9)),1);
    FOI_d = beta_d * Iv / Nh;
    FOI_p = beta_p * (Ip + epsilon_p*Idp) / Nh;
    FOI_v = beta_v * b * (Id + epsilon_dv*Idp) / Nh;
    
    k3 = zeros(12,1);
    k3(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
    k3(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
    k3(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
    k3(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
    k3(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
    k3(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
    k3(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
    k3(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
    k3(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
    k3(10) = Lambda_v - FOI_v*Sv - mu_v*Sv;
    k3(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
    k3(12) = gamma_v*Ev - mu_v*Iv;
    
    y4 = y + dt*k3;
    S=y4(1); Ed=y4(2); Id=y4(3); Rd=y4(4); Ep=y4(5); Ip=y4(6); Rp=y4(7); Idp=y4(8); V=y4(9); Sv=y4(10); Ev=y4(11); Iv=y4(12);
    Nh = max(sum(y4(1:9)),1);
    FOI_d = beta_d * Iv / Nh;
    FOI_p = beta_p * (Ip + epsilon_p*Idp) / Nh;
    FOI_v = beta_v * b * (Id + epsilon_dv*Idp) / Nh;
    
    k4 = zeros(12,1);
    k4(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
    k4(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
    k4(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
    k4(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
    k4(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
    k4(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
    k4(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
    k4(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
    k4(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
    k4(10) = Lambda_v - FOI_v*Sv - mu_v*Sv;
    k4(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
    k4(12) = gamma_v*Ev - mu_v*Iv;
    
    Y_u0(i+1,:) = (y + dt/6*(k1 + 2*k2 + 2*k3 + k4))';
end

%% ==================== SCENARIO 2: u = 0.5 (Constant control) ====================
fprintf('Simulating scenario: u = 0.5...\n');
Y_u01 = zeros(N,12);
Y_u01(1,:) = y0';

for i = 1:N-1
    y = Y_u01(i,:)';
    S=y(1); Ed=y(2); Id=y(3); Rd=y(4); Ep=y(5); Ip=y(6); Rp=y(7); Idp=y(8); V=y(9); Sv=y(10); Ev=y(11); Iv=y(12);
    Nh = max(sum(y(1:9)), 1);
    
    u_const = 0.01;
    bp_eff = beta_p*(1-u_const)*(1-u_const);
    bd_eff = beta_d*(1-u_const)*(1-u_const);
    bv_eff = beta_v*(1-u_const)*(1-u_const);
    Lv_eff = Lambda_v*(1-u_const);
    
    FOI_d = bd_eff * Iv / Nh;
    FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
    FOI_v = bv_eff * b * (Id + epsilon_dv*Idp) / Nh;
    
    k1 = zeros(12,1);
    k1(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
    k1(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
    k1(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
    k1(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
    k1(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
    k1(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
    k1(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
    k1(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
    k1(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
    k1(10) = Lv_eff - FOI_v*Sv - mu_v*Sv;
    k1(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
    k1(12) = gamma_v*Ev - mu_v*Iv;
    
    y2 = y + dt/2*k1;
    S=y2(1); Ed=y2(2); Id=y2(3); Rd=y2(4); Ep=y2(5); Ip=y2(6); Rp=y2(7); Idp=y2(8); V=y2(9); Sv=y2(10); Ev=y2(11); Iv=y2(12);
    Nh = max(sum(y2(1:9)), 1);
    FOI_d = bd_eff * Iv / Nh;
    FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
    FOI_v = bv_eff * b * (Id + epsilon_dv*Idp) / Nh;
    
    k2 = zeros(12,1);
    k2(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
    k2(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
    k2(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
    k2(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
    k2(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
    k2(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
    k2(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
    k2(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
    k2(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
    k2(10) = Lv_eff - FOI_v*Sv - mu_v*Sv;
    k2(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
    k2(12) = gamma_v*Ev - mu_v*Iv;
    
    y3 = y + dt/2*k2;
    S=y3(1); Ed=y3(2); Id=y3(3); Rd=y3(4); Ep=y3(5); Ip=y3(6); Rp=y3(7); Idp=y3(8); V=y3(9); Sv=y3(10); Ev=y3(11); Iv=y3(12);
    Nh = max(sum(y3(1:9)), 1);
    FOI_d = bd_eff * Iv / Nh;
    FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
    FOI_v = bv_eff * b * (Id + epsilon_dv*Idp) / Nh;
    
    k3 = zeros(12,1);
    k3(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
    k3(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
    k3(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
    k3(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
    k3(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
    k3(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
    k3(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
    k3(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
    k3(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
    k3(10) = Lv_eff - FOI_v*Sv - mu_v*Sv;
    k3(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
    k3(12) = gamma_v*Ev - mu_v*Iv;
    
    y4 = y + dt*k3;
    S=y4(1); Ed=y4(2); Id=y4(3); Rd=y4(4); Ep=y4(5); Ip=y4(6); Rp=y4(7); Idp=y4(8); V=y4(9); Sv=y4(10); Ev=y4(11); Iv=y4(12);
    Nh = max(sum(y4(1:9)), 1);
    FOI_d = bd_eff * Iv / Nh;
    FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
    FOI_v = bv_eff * b * (Id + epsilon_dv*Idp) / Nh;
    
    k4 = zeros(12,1);
    k4(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
    k4(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
    k4(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
    k4(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
    k4(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
    k4(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
    k4(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
    k4(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
    k4(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
    k4(10) = Lv_eff - FOI_v*Sv - mu_v*Sv;
    k4(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
    k4(12) = gamma_v*Ev - mu_v*Iv;
    
    Y_u01(i+1,:) = (y + dt/6*(k1 + 2*k2 + 2*k3 + k4))';
end

fprintf('Simulating scenario: u = 0.5...\n');
Y_u03 = zeros(N,12);
Y_u03(1,:) = y0';

for i = 1:N-1
    y = Y_u03(i,:)';
    S=y(1); Ed=y(2); Id=y(3); Rd=y(4); Ep=y(5); Ip=y(6); Rp=y(7); Idp=y(8); V=y(9); Sv=y(10); Ev=y(11); Iv=y(12);
    Nh = max(sum(y(1:9)), 1);
    
    u_const = 0.02;
    bp_eff = beta_p*(1-u_const)*(1-u_const);
    bd_eff = beta_d*(1-u_const)*(1-u_const);
    bv_eff = beta_v*(1-u_const)*(1-u_const);
    Lv_eff = Lambda_v*(1-u_const);
    
    FOI_d = bd_eff * Iv / Nh;
    FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
    FOI_v = bv_eff * b * (Id + epsilon_dv*Idp) / Nh;
    
    k1 = zeros(12,1);
    k1(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
    k1(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
    k1(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
    k1(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
    k1(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
    k1(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
    k1(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
    k1(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
    k1(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
    k1(10) = Lv_eff - FOI_v*Sv - mu_v*Sv;
    k1(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
    k1(12) = gamma_v*Ev - mu_v*Iv;
    
    y2 = y + dt/2*k1;
    S=y2(1); Ed=y2(2); Id=y2(3); Rd=y2(4); Ep=y2(5); Ip=y2(6); Rp=y2(7); Idp=y2(8); V=y2(9); Sv=y2(10); Ev=y2(11); Iv=y2(12);
    Nh = max(sum(y2(1:9)), 1);
    FOI_d = bd_eff * Iv / Nh;
    FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
    FOI_v = bv_eff * b * (Id + epsilon_dv*Idp) / Nh;
    
    k2 = zeros(12,1);
    k2(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
    k2(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
    k2(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
    k2(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
    k2(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
    k2(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
    k2(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
    k2(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
    k2(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
    k2(10) = Lv_eff - FOI_v*Sv - mu_v*Sv;
    k2(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
    k2(12) = gamma_v*Ev - mu_v*Iv;
    
    y3 = y + dt/2*k2;
    S=y3(1); Ed=y3(2); Id=y3(3); Rd=y3(4); Ep=y3(5); Ip=y3(6); Rp=y3(7); Idp=y3(8); V=y3(9); Sv=y3(10); Ev=y3(11); Iv=y3(12);
    Nh = max(sum(y3(1:9)), 1);
    FOI_d = bd_eff * Iv / Nh;
    FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
    FOI_v = bv_eff * b * (Id + epsilon_dv*Idp) / Nh;
    
    k3 = zeros(12,1);
    k3(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
    k3(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
    k3(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
    k3(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
    k3(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
    k3(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
    k3(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
    k3(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
    k3(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
    k3(10) = Lv_eff - FOI_v*Sv - mu_v*Sv;
    k3(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
    k3(12) = gamma_v*Ev - mu_v*Iv;
    
    y4 = y + dt*k3;
    S=y4(1); Ed=y4(2); Id=y4(3); Rd=y4(4); Ep=y4(5); Ip=y4(6); Rp=y4(7); Idp=y4(8); V=y4(9); Sv=y4(10); Ev=y4(11); Iv=y4(12);
    Nh = max(sum(y4(1:9)), 1);
    FOI_d = bd_eff * Iv / Nh;
    FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
    FOI_v = bv_eff * b * (Id + epsilon_dv*Idp) / Nh;
    
    k4 = zeros(12,1);
    k4(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
    k4(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
    k4(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
    k4(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
    k4(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
    k4(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
    k4(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
    k4(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
    k4(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
    k4(10) = Lv_eff - FOI_v*Sv - mu_v*Sv;
    k4(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
    k4(12) = gamma_v*Ev - mu_v*Iv;
    
    Y_u03(i+1,:) = (y + dt/6*(k1 + 2*k2 + 2*k3 + k4))';
end

%% ==================== SCENARIO 3: Optimal control ====================
fprintf('Simulating scenario: Optimal control...\n');

u = 0.5 * ones(N,4);  % u1, u2, u3, u4 initial guess

max_iter = 100;
tol = 1e-4;
rho = 0.5;

for iter = 1:max_iter
    u_old = u;
    
    % Forward solve
    Y_opt = zeros(N,12);
    Y_opt(1,:) = y0';
    
    for i = 1:N-1
        y = Y_opt(i,:)';
        S=y(1); Ed=y(2); Id=y(3); Rd=y(4); Ep=y(5); Ip=y(6); Rp=y(7); Idp=y(8); V=y(9); Sv=y(10); Ev=y(11); Iv=y(12);
        Nh = max(sum(y(1:9)), 1);
        
        u1 = u(i,1); u2 = u(i,2); u3 = u(i,3); u4 = u(i,4);
        bp_eff = beta_p*(1-u1)*(1-u2);
        bd_eff = beta_d*(1-u3)*(1-u4);
        bv_eff = beta_v*(1-u3)*(1-u4);
        Lv_eff = Lambda_v*(1-u4);
        
        FOI_d = bd_eff * Iv / Nh;
        FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
        FOI_v = bv_eff * b * (Id + epsilon_dv*Idp) / Nh;
        
        k1 = zeros(12,1);
        k1(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
        k1(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
        k1(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
        k1(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
        k1(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
        k1(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
        k1(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
        k1(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
        k1(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
        k1(10) = Lv_eff - FOI_v*Sv - mu_v*Sv;
        k1(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
        k1(12) = gamma_v*Ev - mu_v*Iv;
        
        y2 = y + dt/2*k1;
        S=y2(1); Ed=y2(2); Id=y2(3); Rd=y2(4); Ep=y2(5); Ip=y2(6); Rp=y2(7); Idp=y2(8); V=y2(9); Sv=y2(10); Ev=y2(11); Iv=y2(12);
        Nh = max(sum(y2(1:9)), 1);
        FOI_d = bd_eff * Iv / Nh;
        FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
        FOI_v = bv_eff * b * (Id + epsilon_dv*Idp) / Nh;
        
        k2 = zeros(12,1);
        k2(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
        k2(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
        k2(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
        k2(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
        k2(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
        k2(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
        k2(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
        k2(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
        k2(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
        k2(10) = Lv_eff - FOI_v*Sv - mu_v*Sv;
        k2(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
        k2(12) = gamma_v*Ev - mu_v*Iv;
        
        y3 = y + dt/2*k2;
        S=y3(1); Ed=y3(2); Id=y3(3); Rd=y3(4); Ep=y3(5); Ip=y3(6); Rp=y3(7); Idp=y3(8); V=y3(9); Sv=y3(10); Ev=y3(11); Iv=y3(12);
        Nh = max(sum(y3(1:9)), 1);
        FOI_d = bd_eff * Iv / Nh;
        FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
        FOI_v = bv_eff * b * (Id + epsilon_dv*Idp) / Nh;
        
        k3 = zeros(12,1);
        k3(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
        k3(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
        k3(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
        k3(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
        k3(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
        k3(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
        k3(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
        k3(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
        k3(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
        k3(10) = Lv_eff - FOI_v*Sv - mu_v*Sv;
        k3(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
        k3(12) = gamma_v*Ev - mu_v*Iv;
        
        y4 = y + dt*k3;
        S=y4(1); Ed=y4(2); Id=y4(3); Rd=y4(4); Ep=y4(5); Ip=y4(6); Rp=y4(7); Idp=y4(8); V=y4(9); Sv=y4(10); Ev=y4(11); Iv=y4(12);
        Nh = max(sum(y4(1:9)), 1);
        FOI_d = bd_eff * Iv / Nh;
        FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
        FOI_v = bv_eff * b * (Id + epsilon_dv*Idp) / Nh;
        
        k4 = zeros(12,1);
        k4(1) = Lambda_h - FOI_d*S - FOI_p*S - mu_h*S - alpha_d*S - alpha_p*S;
        k4(2) = FOI_d*(S + sigma_p*Rp + sigma_v*V) - (gamma_d + mu_h)*Ed - rho_dp*FOI_p*Ed;
        k4(3) = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id - rho_dp*FOI_p*Id;
        k4(4) = tau_d*Id - mu_h*Rd - phi_d*FOI_p*Rd;
        k4(5) = FOI_p*(S + sigma_d*Rd + sigma_v*V) - (gamma_p + mu_h)*Ep - rho_pd*FOI_d*Ep;
        k4(6) = gamma_p*Ep - (mu_h + delta_p + tau_p)*Ip - rho_pd*FOI_d*Ip;
        k4(7) = tau_p*Ip - mu_h*Rp - phi_p*FOI_d*Rp;
        k4(8) = rho_dp*FOI_p*(Ed+Id) + rho_pd*FOI_d*(Ep+Ip) + phi_d*FOI_p*Rd + phi_p*FOI_d*Rp - (mu_h+delta_dp+tau_dp)*Idp;
        k4(9) = alpha_d*S + alpha_p*S - mu_h*V - sigma_v*FOI_d*V - sigma_v*FOI_p*V;
        k4(10) = Lv_eff - FOI_v*Sv - mu_v*Sv;
        k4(11) = FOI_v*Sv - (gamma_v + mu_v)*Ev;
        k4(12) = gamma_v*Ev - mu_v*Iv;
        
        Y_opt(i+1,:) = (y + dt/6*(k1 + 2*k2 + 2*k3 + k4))';
    end
    
    % Backward solve
    L = zeros(N,12);
    
    for i = N-1:-1:1
        S=Y_opt(i,1); Ed=Y_opt(i,2); Id=Y_opt(i,3); Rd=Y_opt(i,4); Ep=Y_opt(i,5); Ip=Y_opt(i,6);
        Rp=Y_opt(i,7); Idp=Y_opt(i,8); V=Y_opt(i,9); Sv=Y_opt(i,10); Ev=Y_opt(i,11); Iv=Y_opt(i,12);
        Nh = max(sum(Y_opt(i,1:9)), 1);
        
        u1 = u(i,1); u2 = u(i,2); u3 = u(i,3); u4 = u(i,4);
        bp_eff = beta_p*(1-u1)*(1-u2);
        bd_eff = beta_d*(1-u3)*(1-u4);
        
        FOI_d = bd_eff * Iv / Nh;
        FOI_p = bp_eff * (Ip + epsilon_p*Idp) / Nh;
        FOI_v = beta_v*(1-u3)*(1-u4) * b * (Id + epsilon_dv*Idp) / Nh;
        
        lam = L(i+1,:)';
        l1=lam(1); l2=lam(2); l3=lam(3); l4=lam(4); l5=lam(5); l6=lam(6);
        l7=lam(7); l8=lam(8); l9=lam(9); l10=lam(10); l11=lam(11); l12=lam(12);
        
        dlam1 = -(-FOI_d - FOI_p - mu_h - alpha_d - alpha_p)*l1 - FOI_d*l2 - FOI_p*l5 - (alpha_d+alpha_p)*l9;
        dlam2 = -(-(gamma_d+mu_h) - rho_dp*FOI_p)*l2 - gamma_d*l3 - rho_dp*FOI_p*l8;
        dlam3 = -A1 - (-(mu_h+delta_d+tau_d) - rho_dp*FOI_p)*l3 - tau_d*l4 - rho_dp*FOI_p*l8;
        dlam4 = -(-mu_h - phi_d*FOI_p)*l4 - phi_d*FOI_p*l8 - FOI_p*l5;
        dlam5 = -(-(gamma_p+mu_h) - rho_pd*FOI_d)*l5 - gamma_p*l6 - rho_pd*FOI_d*l8;
        dlam6 = -A2 - (-(mu_h+delta_p+tau_p) - rho_pd*FOI_d)*l6 - tau_p*l7 - rho_pd*FOI_d*l8;
        dlam7 = -(-mu_h - phi_p*FOI_d)*l7 - phi_p*FOI_d*l8 - sigma_p*FOI_d*l2;
        dlam8 = -A3 - (-(mu_h+delta_dp+tau_dp))*l8;
        dlam9 = -(-mu_h - sigma_v*FOI_d - sigma_v*FOI_p)*l9 - sigma_v*FOI_d*l2 - sigma_v*FOI_p*l5;
        dlam10 = -(-FOI_v - mu_v)*l10 - FOI_v*l11;
        dlam11 = -(-(gamma_v+mu_v))*l11 - gamma_v*l12;
        dlam12 = -bd_eff*S/Nh*l1 - bd_eff*(S+sigma_p*Rp+sigma_v*V)/Nh*l2 - mu_v*l12;
        
        dlam = [dlam1; dlam2; dlam3; dlam4; dlam5; dlam6; dlam7; dlam8; dlam9; dlam10; dlam11; dlam12];
        L(i,:) = (lam -dt*dlam)';
    end
    
    % Update controls
    for i = 1:N
        S=Y_opt(i,1); Id=Y_opt(i,3); Ip=Y_opt(i,6); Idp=Y_opt(i,8); Sv=Y_opt(i,10);
        Nh = max(sum(Y_opt(i,1:9)), 1);
        
        l1 = L(i,1); l2 = L(i,2); l5 = L(i,5); l10 = L(i,10); l11 = L(i,11);
        
        term_p = (Ip + epsilon_p*Idp) / Nh;
        term_d = Y_opt(i,12) / Nh;
        
        u1_tilde = ((l1 - l2 - l5) * beta_p * (1-u(i,2)) * term_p * S) / B1;
        u2_tilde = ((l1 - l2 - l5) * beta_p * (1-u(i,1)) * term_p * S) / B2;
        u3_tilde = ((l1 - l2 - l5) * beta_d * term_d * S + (l10 - l11) * beta_v * b * (Id + epsilon_dv*Idp)/Nh * Sv) / B3;
        u4_tilde = ((l1 - l2 - l5) * beta_d * term_d * S + (l10 - l11) * beta_v * b * (Id + epsilon_dv*Idp)/Nh * Sv + l10 * Lambda_v) / B4;
        
        u_star = [min(1, max(0, u1_tilde));
                  min(1, max(0, u2_tilde));
                  min(1, max(0, u3_tilde));
                  min(1, max(0, u4_tilde))];
        
        u(i,:) = (1-rho)*u_old(i,:) + rho*u_star';
    end
    
    err = max(abs(u(:) - u_old(:)));
    if err < tol
        fprintf('Optimal control converged after %d iterations!\n', iter);
        break;
    end
end


%% FIGURE: Infected populations AND Cost over time (Two subplots)
figure('Position', [100, 100, 1200, 500]);

% LEFT PLOT: Infected populations
%subplot(1,2,1);
%plot(t, Y_u0(:,3), 'r--', 'LineWidth', 2); hold on;
%plot(t, Y_u01(:,3), 'b-.', 'LineWidth', 2);
%plot(t, Y_u03(:,3), 'g-.', 'LineWidth', 2);
%plot(t, Y_opt(:,3), 'g', 'LineWidth', 2);

%xlabel('Time (days)', 'FontSize', 14);
%ylabel('Infectious population (I_d)', 'FontSize', 14);
%title('Infected Population', 'FontSize', 14);
%legend('u=0(No Control)', 'u=0.01', 'Optimal Control', 'FontSize', 14, 'Location', 'best');
%grid on;
%set(gca, 'FontSize', 14, 'LineWidth', 1.5);

% RIGHT PLOT: Cost over time
%subplot(1,2,2);

% Calculate running cost (integrand)
cost0 = (A1*Y_u0(:,3) + A2*Y_u0(:,6) + A3*Y_u0(:,8));
cost001 = (A1*Y_u01(:,3) + A2*Y_u01(:,6) + A3*Y_u01(:,8) + (B1/2)*0.01^2*4);
cost003 = (A1*Y_u03(:,3) + A2*Y_u03(:,6) + A3*Y_u03(:,8) + (B1/2)*0.03^2*4);

cost_opt = zeros(N,1);
for i = 1:N
    control_cost = (B1/2)*u(i,1)^2 + (B2/2)*u(i,2)^2 + (B3/2)*u(i,3)^2 + (B4/2)*u(i,4)^2;
    cost_opt(i) = (A1*Y_opt(i,3) + A2*Y_opt(i,6) + A3*Y_opt(i,8) + control_cost);
end

% Plot cumulative cost (integral over time)
cum_cost0 = cumsum(cost0 * dt);
cum_cost001 = cumsum(cost001 * dt);
cum_cost003 = cumsum(cost003 * dt);
cum_cost_opt = cumsum(cost_opt * dt);

plot(t, cum_cost0, 'r--', 'LineWidth', 2); hold on;
plot(t, cum_cost001, 'b-.', 'LineWidth', 2);
%plot(t, cum_cost003, 'g-.', 'LineWidth', 2);
plot(t, cum_cost_opt, 'g', 'LineWidth', 2);

xlabel('Time (days)', 'FontSize', 14);
ylabel('Cumulative Cost', 'FontSize', 14);
%title('Associated Cost', 'FontSize', 14);
legend('u=0 (No Control)', 'u=0.01', 'Optimal Control', 'FontSize', 14, 'Location', 'best');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 1.5);


 