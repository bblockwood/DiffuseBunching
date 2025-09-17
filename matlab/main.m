% This script estimates the model on a simulated data set. 

% Set output directory
OUTDIR = '../output/';

% Tax parameters 
k = 300;
t0 = 0.1;
t1 = 0.2;
dT = 0; 

%% ---------Estimation for Bunching around a Kink---------------------------------------------- 
% Import sample income data
simDataBaseline = readtable('../data/sample_kink.csv');

zData = simDataBaseline.income_bin;
hData = table2array(simDataBaseline(:,2));

% Specify starting parameters
mu0 = 5;
e0 = 0.2;
polyDegree = 3; % degree of polynomial for type density

% Estimate best-fit values of e and mu
[hHat, pHat, pSE] = ...
    estimate_model(zData, hData, k, t0, t1, dT, e0, mu0, polyDegree, 1);

% Report estimated parameters
disp(['Estimated elasticity (e): ' num2str(pHat.e)]);
disp(['Estimated lumpiness parameter (mu): ' num2str(pHat.mu)]);

plot(zData, [hData hHat], 'Marker','.','LineWidth',1);
legend('Income histogram', 'Estimated best-fit');
exportgraphics(gcf, fullfile(OUTDIR,'kink_estimation.pdf'));

%% ---------Estimation for Bunching around a Notch---------------------------------------------- 
% Import sample income data
simDataBaseline = readtable('../data/sample_notch.csv');

zData = simDataBaseline.income_bin;
hData = table2array(simDataBaseline(:,2));

% Specify starting parameters
mu0 = 5;
e0 = 0.2;
dT0 = 0.5;
polyDegree = 3; % degree of polynomial for type density

% Estimate best-fit values of e and mu
[hHat, pHat, pSE] = ...
    estimate_model(zData, hData, k, t0, t1, dT0, e0, mu0, polyDegree, 0);

% Report estimated parameters
disp(['Estimated elasticity (e): ' num2str(pHat.e)]);
disp(['Estimated lumpiness parameter (mu): ' num2str(pHat.mu)]);
disp(['Estimated dT: ' num2str(pHat.dT)]);

plot(zData, [hData hHat], 'Marker','.', 'LineWidth',1);
legend('Income histogram', 'Estimated best-fit');
exportgraphics(gcf, fullfile(OUTDIR,'notch_estimation.pdf'));

