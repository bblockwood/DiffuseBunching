function [hHat, thetaFit] = ...
    compute_income_density_quad(z, k, t0, t1, dT, e, mu, theta, hData, hCounterfactual)
% This function plots the model-predicted frequency h at income points z (a vector), when
% there is a bracket threshold at k with marginal tax rates t0 below and t1 above and a
% notch of size dT, for parameters elasticity = e, lumpiness = mu, and with polynomial
% coefficients theta for the income frequency that would be observed under the linear tax
% function T0. 
% 
% This function also has an alternative optional behavior, in which the best-fit
% polynomial coefficients are computed internally using nonlinear least squares, and these
% coefficients are returned in thetaFit. This behavior is triggered by the presence of the
% optional argument hData. In this case only the length of the argument theta is used, to
% determine the order of the polynomial approximation; values of theta are ignored.
% 
% There is also an alternative behavior, triggered by the presence of hCounterfactual, in
% which the density is computed for a specified counterfactual density that would obtain
% under the linear tax function T0. In this case the inputs theta and hData are ignored,
% and hCounterfactual must match z in length. 
% 
% This function works by constructing a 2d matrix with each row representing a different
% income (z) at which we want to compute the observed income density. Within each row,
% columns correspond to different ability levels (n) which contribute to the observed
% income density in that row. Summing across each column produces the total contribution
% to the income density at that income from the different types. We judiciously choose the
% values of the abilities in each column to include key thresholds and to ensure the
% relevant ability range is sufficiently wide to ensure that types beyond that range
% contribute negligible density toward the income in question. 

% Set computational parameters that affect the precision of results
nGridSize = 10000; % number of types to integrate over
gMin = 1e-7;  % definition of "negligible" for choosing n-bounds
z = z(:); % make sure z is a column vector

%% Find income range across which to integrate for distortions
% Compute the income distance at which z +- zRange gives gMin. No point in integrating
% across types more distant than this for distortion at a given z'.
zRange = -(mu/2) * log(gMin*mu);

% convert elasticity to level response: zeta = dz / d(1-t)
zeta = e * k / (1 - t0);
tDiff = t1 - t0;
dZ = zeta * (-tDiff); % z^*_1 - z^*_0 = dZ, negative if tax increases

% In quadratic indirect utility functions, coefficient a on z^2 is same for v0 and v1
a = -1/(2*zeta); % coefficient on z^2 in v0 and v1 are the same

%% Compute distortion at incomes below bracket threshold k

isBelow = (z <= k);
zzLo = z(isBelow); % vector of incomes below bracket threshold k to compute density at

% In this section, index type n by their target income under T0, z^*_0. Ending n for
% each z' is z' + zRange = the n sufficiently high that z^*_0(n) - z' gives gMin.
nLoEnd = zzLo + zRange;

if dT == 0 % pure kink

    % Index type n by target income under T0. At each z', find the n such that v0(z',n) =
    % v0(k,n). At n lower than this, g(z|n) = g0(z|n), so no reason to consider lower n.
    % With quadratic indirect utility, this is just the halfway point between z' and k.
    nLoStart = (zzLo + k)/2;
    nLoStart = min(nLoStart, nLoEnd); % if above nLoEnd, replace with nLoEnd

    % Create grid of n's over which to integrate to get density distortion at each income
    nMatLo = nLoStart + (nLoEnd-nLoStart)./nGridSize * (0:nGridSize);

    Z0UpperBar = nMatLo + (nMatLo - zzLo); % same distance above z^*_0 as z' is below it

    % To compute Z1UpperBar, construct quadratic indirect utility functions v0 and v1 at
    % each z' and n: v0 = a0*z^2 + b0*z + c0, v1 = a1*z^2 + b1*z + c1. 
    b0 = -2*a*nMatLo; % satisfies v'_0(z^*_0) = 0
    b1 = b0 - tDiff; % satisfies v'_1(k) = v'_1(k) - tDiff
    c0 = -a*k^2 - b0*k; % satisfies v_0(k) = 0 (this is a normalization)
    c1 = c0 - dT + tDiff*k;

    v0z = a.*zzLo.^2 + b0.*zzLo + c0; % v0(z')

    % Z1UpperBar is the z s.t. v1(z) - v0(z') = 0, using quadratic formula
    Z1UpperBar = (-b1 - sqrt(b1.^2 - 4*a*(c1 - v0z))) ./ (2*a);

    deltaRLo = Z1UpperBar - Z0UpperBar; % negative, because dominating region shrinks 


elseif dT > 0 % low side of notch is preferable

    % Again, at each z' the lowest n whose choice begins being distorted by the threshold
    % is the halfway point between z' and k.
    nLoStart = (zzLo + k)/2;
    nLoStart = min(nLoStart, nLoEnd); % if above nLoEnd, replace with nLoEnd
    % Create grid of n's over which to integrate to get density distortion at each income
    nMatLo = nLoStart + (nLoEnd - nLoStart) ./ nGridSize * (0:nGridSize);

    % For numerical accuracy, it is useful to add some key "threshold" ability points to
    % this matrix. One such ability threshold is the n s.t. v0(z',n) = v1(k,n). 
    n1ToAdd = zeta ./ (zzLo - k) .* (a*(k^2 - zzLo.^2) - dT);
    
    % For values of z' close to k, this n value will "explode", eventually becoming larger
    % than our already-computed n upper bound. We don't want to add such values to the
    % matrix, so replace them with the previously computed n upper bound value.
    n1ToAdd = max( min(n1ToAdd, nLoEnd), nLoStart); % constrain to existing n interval

    % A second key ability threshold is the ability s.t. v0(z',n) = v1(z^*_1(n),n), since
    % that's the point where a discontinuous dominating income region "disappears" to the
    % right of k. Compute this by solving the resulting quadratic equation in n, with
    % coefficients: 
    an = a + 1/zeta;
    bn = 2*a*dZ + (dZ - zzLo)./zeta - tDiff;
    cn = a*(dZ^2 - zzLo.^2) - dT + tDiff*(k - dZ);

    n2ToAdd1 = (-bn + sqrt(bn.^2 - 4*an*cn)) ./ (2*an); % quadratic formula
    n2ToAdd2 = (-bn - sqrt(bn.^2 - 4*an*cn)) ./ (2*an); 
    n2ToAdd1(n2ToAdd1 < k) = nLoEnd(n2ToAdd1 < k); % if below k, use upper bound instead
    n2ToAdd2(n2ToAdd2 < k) = nLoEnd(n2ToAdd2 < k);
    n2ToAdd1(imag(n2ToAdd1) ~= 0) = nLoEnd(imag(n2ToAdd1) ~= 0); % if imag., use upper bnd
    n2ToAdd2(imag(n2ToAdd2) ~= 0) = nLoEnd(imag(n2ToAdd2) ~= 0);
    n2ToAdd1 = max( min(n2ToAdd1, nLoEnd), nLoStart); % constrain to existing n interval
    n2ToAdd2 = max( min(n2ToAdd2, nLoEnd), nLoStart);

    % Now add these extra n-values to our matrix of types. 
    nMatLo = sort([nMatLo n1ToAdd n2ToAdd1 n2ToAdd2], 2); % sort n values into position
    
    % To compute Z1UpperBar, construct quadratic coefficients for n-matrix
    b0 = -2*a*nMatLo; 
    b1 = b0 - tDiff;
    c0 = -a*k^2 - b0*k;
    c1 = c0 - dT + tDiff*k;

    Z0UpperBar = nMatLo + (nMatLo - zzLo); % same distance above z^*_0 as z' is below it

    v0z = a.*zzLo.^2 + b0.*zzLo + c0;

    Z1UpperBar = (-b1 - sqrt(b1.^2 - 4*a*(c1 - v0z))) ./ (2*a);
    Z1LowerBar = (-b1 + sqrt(b1.^2 - 4*a*(c1 - v0z))) ./ (2*a);
    Z1UpperBar(imag(Z1UpperBar) ~= 0) = nan; % if imaginary, set to nan
    Z1LowerBar(imag(Z1LowerBar) ~= 0) = nan;

    % Compute change in dominating income region, negative when region shrinks
    deltaRLo = max(Z1UpperBar,k) - Z0UpperBar; % if Z1UpperBar is nan, replace with k

    % subtract "hole" if it exists
    hasHoleLo = Z1LowerBar > k;
    deltaRLo(hasHoleLo) = deltaRLo(hasHoleLo) - (Z1LowerBar(hasHoleLo) - k);

elseif dT < 0 % high side of notch is preferable

    % At each z', the minimal n such that the type's density starts being distorted is the
    % n such that v0(z',n) = v1(k,n).
    nLoStart = zeta./(k - zzLo) .* (a*(zzLo.^2 - k^2) + dT);

    % When z' is close to k, this value becomes very negative, so instead use z' - zRange
    % as the starting point.
    nLoStart = max(nLoStart, zzLo - zRange); 
    nLoStart = min(nLoStart, nLoEnd); % if above nLoEnd, replace with nLoEnd

    % Create grid of n's over which to integrate delta R at each (z,n) combination
    nMatLo = nLoStart + (nLoEnd-nLoStart)./nGridSize * (0:nGridSize);

    % Add in the n such that z' = z^*_0(n), which is just n = z' given our index. If
    % that's lower than the n lower bound, replace it with the lower bound.
    nMatLo = sort([nMatLo max(zzLo,nLoStart)], 2); % add these types to n matrix

    % Compute coefficients of quadratic indirect utility functions
    b0 = -2*a*nMatLo;
    b1 = b0 - tDiff;
    c0 = -a*k^2 - b0*k;
    c1 = c0 - dT + tDiff*k;

    % Now compute change in dominating region: Z0UpperBar - Z1UpperBar > 0
    Z0Bar = nMatLo + (nMatLo - zzLo); % same distance on other side of z^*_0
    isHigher = Z0Bar > zzLo;
    Z0UpperBar = Z0Bar;
    zzLoMat = repmat(zzLo, 1, size(Z0Bar,2));
    Z0UpperBar(~isHigher) = zzLoMat(~isHigher);

    v0z = a.*zzLo.^2 + b0.*zzLo + c0; % v0(z')
    % Z1UpperBar is the z s.t. v1(z) - v0(z') = 0, using quadratic formula
    Z1UpperBar = (-b1 - sqrt(b1.^2 - 4*a*(c1 - v0z))) ./ (2*a);

    deltaRLo = Z1UpperBar - Z0UpperBar; % positive if dominating region grows

    % If there is a "hole", subtract that from deltaR
    hasHoleLo = Z0UpperBar < k; % then this type has a hole
    deltaRLo(hasHoleLo) = deltaRLo(hasHoleLo) - (k - Z0UpperBar(hasHoleLo));

end


%% Compute distortion at incomes above bracket threshold k
isAbove = (z > k);
zzHi = z(isAbove); % vector of incomes above bracket threshold k

% In this section we index each type n by their target income under T1.

% In this section, index type n by their target income under T1, z^*_1. Starting n for
% each z' is z' - zRange = the n sufficiently low that z' - z^*_1(n) gives gMin.
nHiStart = zzHi - zRange;

if dT == 0 % pure kink

    % Again find the n halfway between k and z'.
    nHiEnd = (k + zzHi)/2;
    nHiEnd = max(nHiEnd, nHiStart); % if less than nHiStart, replace with nHiStart

    % Create grid of n's over which to integrate delta R at each (z,n) combination
    nMatHi = nHiStart + (nHiEnd-nHiStart)./nGridSize * (0:nGridSize);

    % Now compute change in dominating region: Z1LowerBar - Z0LowerBar < 0
    Z1LowerBar = nMatHi - (zzHi - nMatHi); % same distance below z^*_1 as z' is above it

    % Again compute coefficients of quadratic indirect utility functions
    b1 = -2*a*nMatHi;
    b0 = b1 + tDiff;
    c0 = -a*k^2 - b0*k;
    c1 = c0 - dT + tDiff*k;

    v1Z1 = a.*Z1LowerBar.^2 + b1.*Z1LowerBar + c1; % v1(Z1LowerBar)

    % Z0LowerBar is the z s.t. v0(z) - v1(Z1LowerBar) = 0, using quadratic formula
    Z0LowerBar = (-b0 + sqrt(b0.^2 - 4*a*(c0 - v1Z1))) ./ (2*a);

    deltaRHi = Z1LowerBar - Z0LowerBar; % negative, because dominating region shrinks

elseif dT > 0 % low side of notch is preferable

    % At each z', the maximal n such that the type's density stops being distorted is the
    % n such that v0(k,n) = v1(z',n).
    nHiEnd = zeta./(k - zzHi) .* (a*(zzHi.^2 - k^2) - dT);

    % When z' is close to k, this value explodes, or (when very close) becomes imaginary,
    % so in those cases use the n such that z + zRange is ZLowerBar, and use that if it is
    % a lower n, or if the baseline calculation has an imaginary component.
    nHiEnd = min(nHiEnd, zzHi + zRange); 
    nHiEnd = max(nHiEnd, nHiStart); % if less than nHiStart, replace with nHiStart

    % Create grid of n's over which to integrate delta R at each (z,n) combination
    nMatHi = nHiStart + (nHiEnd-nHiStart)./nGridSize * (0:nGridSize);

    % Add in the n such that z' = z^*_1(n), which is just n = z' given our index. If
    % that's higher than the n upper bound, replace it with the upper bound.
    nMatHi = sort([nMatHi min(zzHi,nHiEnd)], 2); % add these types to n matrix

    % Compute quadratic coefficients of indirect utility functions
    b1 = -2*a*nMatHi;
    b0 = b1 + tDiff;
    c0 = -a*k^2 - b0*k;
    c1 = c0 - dT + tDiff*k;

    % Now compute change in dominating region: Z1LowerBar - Z0LowerBar < 0
    Z1Bar = nMatHi - (zzHi - nMatHi); % same distance on other side of z^*_1 
    zzHiMat = repmat(zzHi, 1, size(Z1Bar,2));
    isLower = Z1Bar < zzHi;
    Z1LowerBar = Z1Bar;
    Z1LowerBar(~isLower) = zzHiMat(~isLower);

    v1z = a.*zzHi.^2 + b1.*zzHi + c1;
    Z0LowerBar = (-b0 + sqrt(b0.^2 - 4*a*(c0 - v1z))) ./ (2*a);

    deltaRHi = Z1LowerBar - Z0LowerBar; % positive if dominating region grows

    % If there is a "hole", subtract that from deltaR
    hasHoleHi = Z1LowerBar > k; % then this type has a hole
    deltaRHi(hasHoleHi) = deltaRHi(hasHoleHi) - (Z1LowerBar(hasHoleHi) - k);

elseif dT < 0 % high side of notch is preferable

    % At each z', the maximal n such that the type's density stops being distorted is the
    % halfway point between k and z.
    nHiEnd = (k + zzHi)/2;
    nHiEnd = max(nHiEnd, nHiStart); % if less than nHiStart, replace with nHiStart

    % Create grid of n's over which to integrate delta R at each (z,n) combination
    nMatHi = nHiStart + (nHiEnd-nHiStart)./nGridSize * (0:nGridSize);

    % For numerical accuracy, it is useful to add some key "threshold" ability points to
    % this matrix. One such ability threshold is the n s.t. v0(k,n) = v1(z',n). 
    n1ToAdd = zeta ./ (zzHi - k) .* (a*(k^2 - zzHi.^2) + dT);

    % For values of z' close to k, this n value will "explode", eventually becoming larger
    % than our already-computed n upper bound. We don't want to add such values to the
    % matrix, so replace them with the previously computed n upper bound value.
    n1ToAdd = max( min(n1ToAdd, nHiEnd), nHiStart); % constrain to existing n interval

    % A second key ability threshold is the ability s.t. v0(z',n) = v1(z^*_1(n),n), since
    % that's the point where a discontinuous dominating income region "appears" to the
    % right of k. Some of these values will actually have z^*_1(n) below 160, in which
    % case nothing is interesting about them, but they also don't hurt anything, so just
    % leave them in. Compute this by solving the resulting quadratic equation in n, with
    % coefficients: 
    an = a + 1/zeta;
    bn = -( 2*a*dZ + (zzHi + dZ)./zeta - tDiff );
    cn = a*(dZ^2 - zzHi.^2) + dT - tDiff*(k + dZ);    

    n2ToAdd1 = (-bn + sqrt(bn.^2 - 4*an*cn)) ./ (2*an); % quadratic formula
    n2ToAdd2 = (-bn - sqrt(bn.^2 - 4*an*cn)) ./ (2*an);
    n2ToAdd1(n2ToAdd1 > k) = nHiEnd(n2ToAdd1 > k); % if above k, use upper bound
    n2ToAdd2(n2ToAdd2 > k) = nHiEnd(n2ToAdd2 > k);
    n2ToAdd1(imag(n2ToAdd1) ~= 0) = nHiEnd(imag(n2ToAdd1) ~= 0); % if imag., use upper bnd
    n2ToAdd2(imag(n2ToAdd2) ~= 0) = nHiEnd(imag(n2ToAdd2) ~= 0);
    n2ToAdd1 = max( min(n2ToAdd1, nHiEnd), nHiStart); % constrain to existing n interval
    n2ToAdd2 = max( min(n2ToAdd2, nHiEnd), nHiStart);

    % Add in the n such that z' = z^*_0(n), which is just n = z' given our index. If
    % that's higher than the n upper bound, replace it with the upper bound.
    nMatHi = sort([nMatHi n1ToAdd n2ToAdd1 n2ToAdd2], 2); % add these types to n matrix

    b1 = -2*a*nMatHi;
    b0 = b1 + tDiff;
    c0 = -a*k^2 - b0*k;
    c1 = c0 - dT + tDiff*k;

    Z1LowerBar = nMatHi - (zzHi - nMatHi); % same distance below z^*_1 as z' is above it

    v1z = a.*zzHi.^2 + b1.*zzHi + c1;

    Z0LowerBar = (-b0 + sqrt(b0.^2 - 4*a*(c0 - v1z))) ./ (2*a);
    Z0UpperBar = (-b0 - sqrt(b0.^2 - 4*a*(c0 - v1z))) ./ (2*a);
    Z0LowerBar(imag(Z0LowerBar) ~= 0) = nan;
    Z0UpperBar(imag(Z0UpperBar) ~= 0) = nan;

    % Compute change in dominating income region, negative when region shrinks
    deltaRHi = Z1LowerBar - min(Z0LowerBar,k); % if Z0LowerBar is nan, replace with k

    % subtract "hole" if it exists
    hasHoleHi = Z0UpperBar < k;
    deltaRHi(hasHoleHi) = deltaRHi(hasHoleHi) - (k - Z0UpperBar(hasHoleHi));

end


% Compute bunching distortion
lambda = 1/mu;

% Compute distortion at each z' and n relative to h0 and h1 under linear tax
bMatLo = lambda*exp(-2*lambda*abs(nMatLo-zzLo)).*(exp(-lambda*deltaRLo)-1);
bMatHi = lambda*exp(-2*lambda*abs(zzHi-nMatHi)).*(exp(-lambda*deltaRHi)-1);
% Create matrices of types for computing density under linear tax T0
n0MatLo = nMatLo;
n0MatHi = nMatHi-dZ;


%% Compute income density under manually specified counterfactual

if exist('hCounterfactual', 'var') && ~isempty(hCounterfactual)
    hCounterfactual = hCounterfactual(:);
    assert(length(hCounterfactual) == length(z), "hCounterfactual must be the same length as z");

    % use interpolation to find counterfactual type density f(n) at each n in matrices
    % TODO: this should really be computed using the density of *types*, not the density over observed incomes -- could use FFT for this. 
    f0MatLo = interp1(z, hCounterfactual, n0MatLo, 'linear', 'extrap');
    f0MatHi = interp1(z, hCounterfactual, n0MatHi, 'linear', 'extrap');

    % compute distortion at each z' and n relative to h0 and h1 under linear tax
    bf0MatLo = bMatLo.*f0MatLo;
    bf0MatHi = bMatHi.*f0MatHi;

    % integrate across each row of hMatLo
    bLo = nan(size(bf0MatLo, 1), 1);
    bHi = nan(size(bf0MatHi, 1), 1);
    for row = 1:size(bf0MatLo, 1)
        bLo(row) = trapz(n0MatLo(row,:), bf0MatLo(row,:));
    end
    for row = 1:size(bf0MatHi, 1)
        bHi(row) = trapz(n0MatHi(row,:), bf0MatHi(row,:));
    end

    % compute density under counterfactual
    h0CFLo = hCounterfactual(isBelow);
    h1CFHi = interp1(z-dZ, hCounterfactual, zzHi, 'linear', 'extrap');
    hHat = [bLo; bHi] + [h0CFLo; h1CFHi];

    thetaFit = [];
    return;
    
end

%% Compute income density under polynomial counterfactual and estimate theta coefficients

% We use the fact that h(z) is the integral g(z|n)*f(n) over n, where
% f(n) is a polynomial: f(n) = thetaStar_0 + thetaStar_1*n + thetaStar_2*n^2 + ... The
% integral is computed using the trapezoidal method: h(z) ~=~ sum ((n_{k+1} - n_k) *
% (g(z|n_k)*f(n_k) + g(z|n_{k+1})*f(n_{k+1}))/2 ). By writing out the polynomial form of
% f(n), this can be factored out to give a linear expression in the coefficients
% thetaStar_0, thetaStar_1, thetaStar_2, ..., which means that, if desired, they can be
% estimated very quickly using least squares minimization with Matlab's backslash
% function. The tricky part is just to compute the vector that is linearly weighted by
% each thetaStar coefficient. That is what is done below, with each Y referring to the
% coefficients.


polyDegree = length(theta) - 1;

nMatLoDiff = nMatLo(:,2:end) - nMatLo(:,1:end-1);
nMatHiDiff = nMatHi(:,2:end) - nMatHi(:,1:end-1);

% Loop over polynomial degree to construct the weight on each coefficient
YMat = nan(length(z),polyDegree+1);
for p = 0:polyDegree
    bnpLo = ( bMatLo(:,1:end-1).*n0MatLo(:,1:end-1).^p + ... 
        bMatLo(:,2:end).*n0MatLo(:,2:end).^p )./2;
    bnpHi = ( bMatHi(:,1:end-1).*n0MatHi(:,1:end-1).^p + ...
        bMatHi(:,2:end).*n0MatHi(:,2:end).^p )./2;

    YpLo = sum(nMatLoDiff.*bnpLo, 2) + zzLo.^p;
    YpHi = sum(nMatHiDiff.*bnpHi, 2) + (zzHi-dZ).^p;
    % Add further adjustment to account for star vs. not

    Yp = [YpLo; YpHi]; % construct full vector below and above k
    YMat(:,p+1) = Yp; % add vector
end

% Estimate best-fit theta polynomial coefficients
if exist('hData', 'var') && ~isempty(hData)
    assert(length(hData) == length(z), "hData must be the same length as z");

    % Solve in a more stable numerical space
    YMat_scaled = YMat / norm(YMat, 'fro'); 
    hData_scaled = hData / norm(hData); 
    thetaFit = ( YMat_scaled \ hData_scaled )'; % Solve in a more stable numerical space
    thetaFit = thetaFit * (norm(hData) / norm(YMat, 'fro')); % Rescale back 

    % adjust theta0 (intercept) so that predicted pdf sum matches data
    thIntAdj = ( sum(hData) - sum(YMat*thetaFit') ) / sum(YMat(:,1));
    thetaFit(1) = thetaFit(1) + thIntAdj;
else 
    thetaFit = theta(:)';
end

% Compute predicted income density
hHat = YMat*thetaFit';
