clc; clear; close all;

%% 1. Dengue Data (365 days)
Id_data = [52 20 23 29 18 8 22 26 7 41 16 18 10 35 49 42 31 38 33 9 36 24 35 19 29 14 25 43 53 49 60 42 45 22 47 14 82 48 83 65 23 79 86 84 101 67 46 17 90 115 78 70 120 62 21 22 147 94 114 83 62 317 26 19 28 38 288 108 159 169 249 234 244 212 248 151 352 329 392 394 326 195 159 262 383 429 386 416 358 204 294 317 492 425 406 337 138 391 420 330 375 321 375 114 394 429 364 444 319 280 164 331 409 394 393 386 278 138 209 343 395 319 428 408 190 325 448 434 382 325 252 134 202 466 380 357 356 311 173 247 430 412 470 430 432 165 367 568 552 473 445 363 158 364 580 573 487 625 586 245 279 685 636 665 622 350 359 382 740 678 664 668 632 219 514 854 735 556 490 396 263 374 1042 782 715 700 781 308 413 953 857 841 758 755 510 619 950 942 814 762 803 468 659 1143 983 1041 964 928 506 651 1162 1147 1101 1069 1034 488 834 1195 1179 912 500 833 460 792 1139 1007 920 788 745 436 593 500 705 633 615 567 500 572 636 610 565 490 565 200 300 516 455 421 377 411 200 397 387 288]';
tspan = (0:1:length(Id_data)-1)';

%% 2. Settings & Initial Conditions
% We use a large population to support a high peak (~1200)
Nh_total = 80000; 
y0_base = [Nh_total*0.95, 0, 0, 0, Nh_total*0.051, 30000, 100, 20]; 

% p = [beta_d   beta_v  alpha_d sigma_v Lambda_h Lambda_v b Ed0 Id0]
% gamma_d/v are set low to delay the peak. Ed0/Id0 are low to slow the start.

params0 = [.6, 01.16, 0.002, 0.02, 20, 10, 100, 30];

lb = [0.58, 01.14, 0.001, 0.01, 10, 1, 40, 5];
ub = [.620, 1.2, 0.003000, .060, 150, 80, 500, 150];

%% 3. Optimization
% Weighted cost: Focus 10x more on the end of the year and 5x on the peak height
options = optimset('Display','iter','MaxIter',3000,'TolX',1e-10,'TolFun',1e-10);
costfun = @(p) weighted_cost(p, y0_base, tspan, Id_data);

% Use fminsearchbnd (Requires the function file in your folder)
params_est = fminsearchbnd(costfun, params0, lb, ub, options);

disp('Estimated Parameters:');
disp(vpa(params_est,9))

%% 4. Final Simulation & Plot
[t, Y] = run_sim(params_est, y0_base, tspan);
Id_model = Y(:,3);

% figure('Color', 'w', 'Name', 'Dengue Model Fit');
plot(tspan, Id_data, 'r', 'LineWidth', 2.0); hold on;
plot(t, Id_model, 'b-', 'LineWidth', 2.5);
xlabel('Days','Interpreter','latex','FontSize',20); ylabel('Infected Individuals ($I_d$)','Interpreter','latex','FontSize',20);
% title('Dengue Fitting: Delayed Peak & Low Tail');
legend({'Real Data','Model Output'},'Location', 'best'); grid on;

%% --- Helper Functions ---

function cost = weighted_cost(p, y0_base, tspan, data)
    [~, Y] = run_sim(p, y0_base, tspan);
    if length(Y(:,3)) < length(tspan), cost = 1e15; return; end
    
    model_vals = Y(:,3);
    weights = ones(size(data));
    
    % Target Peak Area (around day 200-300)
    [~, peakIdx] = max(data);
    weights(peakIdx-30:peakIdx+30) = 5; 
    
    % Target Tail Area (last 50 days) - Force the drop
    weights(end-50:end) = 15; 
    
    cost = sum(weights .* (model_vals - data).^2);
    
    errors = data - model_vals;       % 1. Residuals
squared_errors = errors.^2;        % 2. Square them
mse = mean(squared_errors);        % 3. Mean (MSE)
% rmse_val = sqrt(mse);              % 4. Root (RMSE)
rmse_val = sqrt(mean(errors.^2));

nrmse = rmse_val / (max(data) - min(data))

nrmse_percent = 100 * nrmse;
% fprintf('The RMSE of the model is: %.4f\n', rmse_val);
end

function [t, Y] = run_sim(p, y0_base, tspan)
    y0 = y0_base;
    y0(2) = p(5); % Optimized Ed0
    y0(3) = p(6); % Optimized Id0
    y0(7) = p(7); % Optimized Ev0
    y0(8) = p(8); % Optimized Iv0
    [t, Y] = ode45(@(t,y) dengue_ode(t,y,p), tspan, y0);
end

function dydt = dengue_ode(~,y,p)
    S=y(1); Ed=y(2); Id=y(3); Rd=y(4); V=y(5); Sv=y(6); Ev=y(7); Iv=y(8);
    Nh = S + Ed + Id + Rd + V;
    
    % Parameters
    beta_d=p(1); beta_v=p(2); gamma_d=.15; gamma_v=.11; 
    tau_d=.1428; alpha_d=p(3); sigma_v=p(4); Lambda_h=2.5; 
    Lambda_v=5000; b=.3;
    
    % Hard-coded constants to ensure rapid decline
    delta_d=0.00095; mu_h=0.000038; 
%     mu_v=0.14; % Higher mosquito mortality = faster crash at the end
    mu_v=0.14;
    
    dS  = Lambda_h - beta_d*(Iv/Nh)*S - (mu_h + alpha_d)*S;
    dEd = beta_d*(Iv/Nh)*(S + sigma_v*V) - (gamma_d + mu_h)*Ed;
    dId = gamma_d*Ed - (mu_h + delta_d + tau_d)*Id;
    dRd = tau_d*Id - mu_h*Rd;
    dV  = alpha_d*S - mu_h*V - sigma_v*beta_d*(Iv/Nh)*V;
    
    dSv = Lambda_v - beta_v*b*(Id/Nh)*Sv - mu_v*Sv;
    dEv = beta_v*b*(Id/Nh)*Sv - (gamma_v + mu_v)*Ev;
    dIv = gamma_v*Ev - mu_v*Iv;
    
    dydt = [dS; dEd; dId; dRd; dV; dSv; dEv; dIv];
end