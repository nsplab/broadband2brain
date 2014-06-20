% Cadzow denoising algorithm
% From Blu "Sparse Sampling of Signal Innovations"
% Input and output is vector of samples y

function y_new = cadzow_denoise(y, K)

% Make Toeplitz matrix approx square
N = length(y);
L1 = ceil((N+1)/2);
L2 = floor((N+1)/2);

% Form matrix A
A = zeros(L1, L2);
for i = 1 : L1,
    for j = 1 : L2,
        A(i, j) = y(L2 + i - j);
    end
end

diags = zeros(size(y));

count = 0;
%%{
if min(L1,L2) > K
    
    while 1
    %for i = 1:3  % TESTING: 3 iterations
        
        count = count + 1;
        
        % Reduce dimension using SVD
        [U S V] = svd(A);
        S_prime = zeros(size(S));
        S_prime(1:K,1:K) = S(1:K,1:K);
        A = U*S_prime*V';

        % Average diagonals to produce Toeplitz matrix
        for d = -(L1-1) : L2-1
            val = mean(diag(A,d));
            diags(d+L1) = val;
            for i = max(1, 1-d) : min(L2-d, L1)
                A(i,i+d) = val;
            end
        end
        
        % original terminating condition: 0.0001 (used for simulation)
        % 0.1 gives approx 3 iterations
        if S(K+1,K+1)/S(K,K) < 0.0001  % change terminating condition here
            break;
        end
        
    end

end
%%}

count

y_new = diags(end:-1:1);

end