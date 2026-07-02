%% Multidisease Model Sensitivity Analysis with PRCC
% Self-contained MATLAB code: PRCC, scatterplots, bar chart, heatmap, and time-series

clc; clear; close all;

% Step 1: Define Parameters and Sampling
numSamples = 500; % Latin hypercube samples
paramNames = {'\beta_d','\beta_p','\beta_v','b','\gamma_d','\gamma_p','\tau_d','\tau_p','\mu_h','\gamma_v','\mu_v'};
numParams = length(paramNames);

% Define parameter ranges (example)
lb = [0.1,0.01, 0.05, 0.02, 0.2, 0.1, 0.1, 0.05, 0.05, 0.01 ,0.01];  % lower bounds
ub = [0.9, 0.05,0.5, 0.2, 1.0, 0.8, 0.9, 0.4, 0.6, 0.05,0.05];     % upper bounds

% Latin Hypercube Sampling
lhs = lhsdesign(numSamples, numParams);
params = lhs * (ub - lb) + lb;

%% Step 2: Define Model Output (Example: R0d, R0p, R0v)
% Simplified reproduction number formulas
beta_d = params(:,1); beta_p = params(:,2); beta_v = params(:,3);
b = params(:,4); gamma_d = params(:,5); gamma_p = params(:,6);
tau_d = params(:,7); tau_p = params(:,8); mu_h = params(:,9);gamma_v = params(:,10);mu_v = params(:,11);

R0d = sqrt((beta_d .* beta_v .* b) ./ (mu_h .* (gamma_d + mu_h)));
R0p = (beta_p) ./ ((gamma_p + mu_h) .* (mu_h + tau_p));
R0v = R0d .* 0.7 + R0p .* 0.3; % synthetic combined measure

Outputs = [R0d, R0p, R0v];
outputNames = {'R_{0d}','R_{0p}','R_{0v}'};

%% Step 3: PRCC Calculation
PRCC = zeros(numParams,3);
for j = 1:3
    for i = 1:numParams
        PRCC(i,j) = partialcorr(params(:,i), Outputs(:,j), params(:,setdiff(1:numParams,i)), 'Type','Spearman');
    end
end

%% Step 4: Scatterplots (for first output)
figure('Name','Scatterplots','Position',[100 100 800 600]);
for i = 1:numParams
    subplot(3,3,i);
    scatter(params(:,i), R0d, 30, 'filled', 'MarkerFaceAlpha', 0.5);
    xlabel(paramNames{i}, 'Interpreter','tex');
    ylabel('R_{0d}','Interpreter','tex');
    title(['Scatter: ', paramNames{i}],'FontWeight','normal');
    grid on;
end
sgtitle('Parameter Scatterplots for R_{0d}','FontWeight','bold');

%% Step 5: Vertical Bar Chart (Blue = Negative, Red = Positive)
figure('Name','PRCC Bar Chart','Position',[100 100 800 600]);
for j = 1:3
    subplot(1,3,j);
    posVals = PRCC(:,j) > 0;
    negVals = PRCC(:,j) <= 0;
    hold on;
    bar(find(posVals), PRCC(posVals,j), 'FaceColor',[0.85 0.1 0.1]); % red
    bar(find(negVals), PRCC(negVals,j), 'FaceColor',[0.1 0.3 0.8]); % blue
    hold off;
    set(gca,'XTick',1:numParams,'XTickLabel',paramNames,'XTickLabelRotation',45);
    ylabel('PRCC Value');
    title(['PRCC for ', outputNames{j}]);
    ylim([-1 1]);
    grid on;
end
sgtitle('Partial Rank Correlation Coefficients (PRCC)','FontWeight','bold');

%% Step 6: Heatmap of PRCC
figure('Name','PRCC Heatmap','Position',[100 100 700 500]);
h = heatmap(outputNames, paramNames, PRCC, 'Colormap',redbluecmap);
h.Title = 'Sensitivity Heatmap (PRCC)';
h.XLabel = 'Model Outputs';
h.YLabel = 'Parameters';
h.ColorLimits = [-1 1];

%% Step 7: Time-Series Simulation (synthetic example)
t = linspace(0,100,200);
I_d = 0.05*exp(0.05*t)./(1+0.05*exp(0.05*t));
I_p = 0.03*exp(0.04*t)./(1+0.04*exp(0.04*t));
I_dp = 0.01*exp(0.06*t)./(1+0.06*exp(0.06*t));

figure('Name','Time Series','Position',[100 100 800 500]);
plot(t, I_d, 'r-', 'LineWidth', 2); hold on;
plot(t, I_p, 'b--', 'LineWidth', 2);
plot(t, I_dp, 'k-.', 'LineWidth', 2);
xlabel('Time (days)','FontSize',13);
ylabel('Proportion Infected','FontSize',13);
legend({'Dengue','Pneumonia','Co-infection'}, 'Location','best','FontSize',13);
title('Simulated Time Series of Infections','FontSize',13);
grid on;
