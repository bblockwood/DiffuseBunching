function [hHat, pHat, pSE] = estimate_model(zData, hData, k, t0, t1, ...
    dT0, e0, mu0, polyDegree, dTFixedFlag, hCounterfactual)
% This function estimates the parameters of the model to fit the provided data using
% maximum likelihood. Estimates of e, mu, and (if desired) dT are returned in
% the structure pHat, along with their standard errors in the structure pSE. If the
% optional argument dTFixedFlag is nonzero, then the model estimates only the elasticity,
% not dT. 

% If hCounterfactual is provided, it is used as the counterfactual income density. 

% If dTFixedFlag is not provided, set it to zero (i.e., estimate best-fit dT too)
if ~exist('dTFixedFlag', 'var') 
    dTFixedFlag = 0; 
end
if ~exist('hCounterfactual', 'var') 
    hCounterfactual = []; 
end

% Starting values
theta0 = zeros(1, polyDegree+1); % used only for size of theta
x0 = [e0 mu0 dT0];
lb = [0 0 -10]; ub = [5 mu0*5 10]; % third elements ignored if dTFixedFlag is nonzero

% eStep = 0.001;
eStep = 0.0005*e0;
% muStep = 0.01;
muStep = 0.0005*mu0;
dTStep = 0.01;
stepsize = [eStep muStep dTStep];

function L = compute_log_likelihood(e, mu, dT)
    if exist('hCounterfactual', 'var') && ~isempty(hCounterfactual)
        hHat = compute_income_density(zData, k, t0, t1, dT, e, mu, theta0, hData, hCounterfactual);
        hHat = max(hHat, 1e-10);   
    else
        hHat = compute_income_density(zData, k, t0, t1, dT, e, mu, theta0, hData);
    end
        L = sum(hData .* log(hHat));
end

if dTFixedFlag == 0 % estimate best-fit dT, e, and mu
    obj = @(x) -compute_log_likelihood(x(1), x(2), x(3));
else % estimate only e and mu, with dT fixed at dT0
    x0 = x0(1:2);
    lb = lb(1:2); ub = ub(1:2);
    stepsize = stepsize(1:2);
    obj = @(x) -compute_log_likelihood(x(1), x(2), dT0);
end

options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'interior-point', ...
    'StepTolerance', 1e-100, 'MaxFunEvals', 1e4, ...
    'FiniteDifferenceStepSize', stepsize, 'FiniteDifferenceType', 'central');
problem = createOptimProblem('fmincon', 'objective', obj, 'x0', x0, ...
    'lb', lb, 'ub', ub, 'options', options);
ms = MultiStart('UseParallel', true);
[xSol, fval] = fmincon(problem.objective, ...
    problem.x0, [], [], [], [], problem.lb, problem.ub, [], problem.options);

pHat.fval = fval; % store for comparison, diagnostic
pHat.e = xSol(1);
pHat.mu = xSol(2);

if dTFixedFlag == 0
    pHat.dT = xSol(3);
    dT = pHat.dT;
else
    dT = dT0;
end

if exist('hCounterfactual', 'var') && ~isempty(hCounterfactual)
    [hHat, pHat.theta] = compute_income_density(zData, k, t0, t1, ...
        dT, pHat.e, pHat.mu, theta0, hData, hCounterfactual);
else
    [hHat, pHat.theta] = compute_income_density(zData, k, t0, t1, ...
        dT, pHat.e, pHat.mu, theta0, hData);
end


% Compute covariance matrix for parameter estimates, if they are requested.
if nargout > 2
    % Package negloglikelihood function in the way 'mlecov'expects
    negloglf = @(x,data,cens,freq) obj(x);

    if dTFixedFlag == 0
        pVec = [pHat.e pHat.mu pHat.dT];
    else
        pVec = [pHat.e pHat.mu];
    end
    
    pCov = mlecov(pVec, zData,'nloglf', negloglf, 'Frequency', hData, 'Options', ...
        statset('DerivStep', stepsize));
    
    seVec = sqrt(diag(pCov))';
    pSE.e  = seVec(1);
    pSE.mu = seVec(2);
    if dTFixedFlag == 0
        pSE.dT = seVec(3);
    end
end

end
