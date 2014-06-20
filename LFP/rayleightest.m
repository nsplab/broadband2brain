function [pval,R] = rayleightest(theta,rho)
% function [pval,R] = rayleightest(theta,rho)
%   Rayleigh's test of a uniform circular distribution of unit vectors
%   against a unimodal alternative (Fisher 1993, p. 70). A second input can
%   be used to assign weights to the vectors. The test then is based on the
%   modification of the Rayleigh test by Moore (1980).
%   Theta should be in radians.
if nargin == 1
    n  = length(theta);
    C = sum(cos(theta));
    S = sum(sin(theta));
    R = sqrt(C^2 + S^2)/n;
    Z = n*R^2;
    pval = exp(-Z)*(1+(2*Z - Z^2)/(4*n) - (24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4)/(228*n^2));
else
    n = length(theta);
    theta = reshape(theta,n,1);
    rho = reshape(rho,n,1);
    vects = sortrows([theta rho],2);
    vects(:,2) = [1:n]';
    C = sum(vects(:,2).*cos(vects(:,1)));
    S = sum(vects(:,2).*sin(vects(:,1)));
    R = sqrt(C^2 + S^2)/(n^(3/2));
    psigtest = [0.999 0.500 0.050 0.025 0.010 0.005 0.001];
    [grid_pval,grid_n] = meshgrid(psigtest,[2:10 12:2:30 40:20:100 max(n,10000)]);
    crit_R = [0.354 0.791 1.058 1.060 1.061 1.061 1.061;
              0.014 0.693 1.095 1.124 1.143 1.149 1.154;
              0.020 0.620 1.090 1.146 1.192 1.212 1.238;
              0.023 0.588 1.084 1.152 1.216 1.250 1.298;
              0.022 0.568 1.074 1.152 1.230 1.275 1.345;
              0.021 0.556 1.066 1.150 1.238 1.291 1.373;
              0.020 0.546 1.059 1.148 1.242 1.300 1.397;
              0.020 0.538 1.053 1.146 1.245 1.307 1.416;
              0.020 0.532 1.048 1.144 1.248 1.313 1.432;
              0.019 0.523 1.042 1.140 1.252 1.322 1.456;
              0.019 0.518 1.037 1.136 1.252 1.325 1.470;
              0.019 0.514 1.031 1.132 1.250 1.327 1.480;
              0.019 0.510 1.027 1.129 1.248 1.328 1.487;
              0.019 0.507 1.024 1.127 1.247 1.329 1.492;
              0.019 0.505 1.022 1.126 1.246 1.330 1.496;
              0.019 0.503 1.021 1.125 1.246 1.331 1.499;
              0.019 0.502 1.019 1.124 1.246 1.332 1.501;
              0.019 0.500 1.018 1.124 1.246 1.333 1.502;
              0.018 0.499 1.016 1.123 1.245 1.334 1.502;
              0.018 0.494 1.012 1.119 1.243 1.332 1.504;
              0.018 0.489 1.007 1.115 1.241 1.329 1.506;
              0.018 0.487 1.005 1.113 1.240 1.329 1.508;
              0.018 0.485 1.004 1.112 1.240 1.329 1.509;
              0.018 0.481 0.999 1.109 1.239 1.329 1.517];
    R_crits = interp2(grid_pval,grid_n,crit_R,psigtest,n*ones(1,length(psigtest)),'linear');
    pval = interp1(R_crits,psigtest,R,'linear',NaN);
    if isnan(pval) & (R < min(R_crits)), pval = psigtest(1);
    elseif isnan(pval) & (R > max(R_crits)), pval = psigtest(end); end
end