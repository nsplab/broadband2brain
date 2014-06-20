% Computes E[sigmae | y] by assuming t is correct and integrating out c

function sigmae = sigmae_unbiased(N, K, y, G)
% N = num samples
% K = num spikes
% y = samples
% G(n, k) = h(nT - t_k)

A = sum(y.^2);
Ak = -2*G'*y;
Akk = G'*G;

for k = 1 : K
    % Integrate out c_k
    if Akk(k,k) == 0
        continue;
    end
    A = A - Ak(k)^2/(4*Akk(k,k));
    Ak = Ak - Ak(k)*Akk(:,k)/Akk(k,k);
    Akk = Akk - Akk(:,k)*Akk(k,:)/Akk(k,k);
end

sigmae = sqrt(A/2) * gamma(0.5*(N-K-1)) / gamma(0.5*(N-K));

if ~isreal(sigmae)
    disp('WARNING: sigma_e non-real');
    sigmae = real(sigmae);
end

if isnan(sigmae)
    gamma(0.5*(N-K))
    disp('ERROR: sigma_e = NaN');
end

end