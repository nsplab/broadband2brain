function [t, c, sigma, elapsed_time, t_inter, c_inter, sig_inter gaPhase1_iterTime gaPhase2_iterTime] = reconstruct_GA(y, T, args, prev_data)
% Genetic algorithm gibbs sampler

tic;

sigmah = args.sigmah;
N = args.L;
K = args.K;
I_e = args.I_e;
I_m = args.I_m;
se = args.sigmae;

gaPhase1_iterTime = NaN;
gaPhase2_iterTime = NaN;

y = y';
scalefactor = 1;

t0 = 0; %min(time);  % Take min(time) to be 0
tf = T; %max(time) - min(time);
%dt = get_dt(); %tf / (length(time) - 1);
%time = time - t0;
n = linspace(0, tf, N);  % Vector of nT values for n = 0, 1 ... N-1 (points in time at which to take samples y)

% Initial conditions
c=zeros(1,K);
t=zeros(1,K);

% TESTING: start with correct answer
%c = args.c;
%t = args.t;

% Discrete values of t used to draw from distribution
% 100 is number of values
t_sampled=linspace(0, tf, 100);  % Edited by Alex Wein, changed 20 to tf

% Save intermediate results
c_inter = [];
t_inter = [];
sig_inter = [];

c_best = c;
t_best = t;
f_best = f(t_best, c_best, sigmah, n, y);

% PHASE I

for iter = 1 : I_e
    
    phase1ID = tic;

    % Update c_k

    for k = 1:K

        alpha_0=(1/(2*se^2))*sum(exp(-((n-t(k)).^2)./(sigmah^2)));

        d=zeros(1,length(y));

        for kprime=1:K

            if(kprime ~= k)

                d=d+(c(kprime))*exp(-((n-t(kprime)).^2)./(2.*sigmah^2));

            end

        end



        beta_0=(1/(se^2))*sum((exp(-((n-t(k)).^2)./(2.*sigmah^2))).*(d-y));



        c(k)=normrnd(-beta_0/(2*alpha_0), sqrt(1/(2*alpha_0)));

        % TESTING: set c(k) to mean of distribution
        %c(k) = -beta_0/(2*alpha_0);



    end

    c_inter = [c_inter c'];  % save intermediate results



    % Update t_k

    for k = 1:K



        gamma_0=c(k)^2;  % gamma_k



        d=zeros(1,length(y));

        for kprime=1:K

            if(kprime ~= k)

                d=(c(kprime))*exp(-((n-t(kprime)).^2)./(2.*sigmah^2))+d;

            end

        end



        delta=2*c(k)*(d-y);  % delta is now equal to nu_k in paper (eq 25), note that it takes on a different value for each n

        repn=repmat(n, length(t_sampled), 1);

        repdelta=repmat(delta, length(t_sampled),1);

        rept_sampled=repmat(t_sampled', 1,N);



        lnprobs_t=-1/(2*se^2)*sum(gamma_0*exp(-((repn-rept_sampled).^2)./(sigmah^2))+ ...
            repdelta.*exp(-((repn-rept_sampled).^2)./(2*sigmah^2)),2);  % log probability of each t value in t_sampled

        probs_t=exp(scalefactor*lnprobs_t)/sum(exp(scalefactor*lnprobs_t));  % probability of each t value in t_sampled

        if(norm(probs_t)==0)

            t(k)=rand(1)*tf;  % Edited by Alex Wein, changed 20 to tf
            disp('Warning: t_k chosen randomly');

        else

            t_sampled_nonzero=t_sampled;

            probs_t_nonzeros=probs_t;

            t(k)=discrete_sample(1, t_sampled_nonzero, probs_t_nonzeros);  % sample t_k from discrete values

        end
        
    end

    %t_inter = [t_inter t'];  % save intermediate results

    % Choose next generation
    
    f_new = f(t, c, sigmah, n, y);
    if f_new > f_best
        t_best = t;
        c_best = c;
        f_best = f_new;
    else
        t = t_best;
        c = c_best;
    end
    
    gaPhase1_iterTime = toc(phase1ID);
    
end

% PHASE II

for i = 1 : I_m
    
    phase2ID = tic;
    
    for j = 1 : K
        
        t(j) = unifrnd(0, tf);
        
        % special step for first iteration
        if i == 1
            [val ind] = max(y-calc_z(t, c, sigmah, n));
            t(j) = n(ind);
        end
        
        % Update c_k

        for k = 1:K

            alpha_0=(1/(2*se^2))*sum(exp(-((n-t(k)).^2)./(sigmah^2)));

            d=zeros(1,length(y));

            for kprime=1:K

                if(kprime ~= k)

                    d=d+(c(kprime))*exp(-((n-t(kprime)).^2)./(2.*sigmah^2));

                end

            end



            beta_0=(1/(se^2))*sum((exp(-((n-t(k)).^2)./(2.*sigmah^2))).*(d-y));



            c(k)=normrnd(-beta_0/(2*alpha_0), sqrt(1/(2*alpha_0)));

            % TESTING: set c(k) to mean of distribution
            %c(k) = -beta_0/(2*alpha_0);



        end

        %c_inter = [c_inter c'];  % save intermediate results
        
        f_new = f(t, c, sigmah, n, y);
        if f_new > f_best
            t_best = t;
            c_best = c;
            f_best = f_new;
        else
            t = t_best;
            c = c_best;
        end
        
    end
    
    gaPhase2_iterTime = toc(phase2ID);
    
end


%mtheta=mean(thetaest(end-mean_iter:end,:));

% Original outputs
ckhat=c;
tkhat=t;
%sigmahat=mtheta(end);

% New outputs
c_hat_orig = c;  % c_hat is improved below
t_hat = t + t0;
time_y = n + t0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hmatrix=zeros(N,K);

for k = 1:K

    for nn=1:N

        Hmatrix(nn,k)=gaussian(tkhat(k)-n(nn), args);

    end

end

ckhat2=Hmatrix\y';



% Output improved c_hat
if all(isnan(ckhat2)==0)
    c_hat = ckhat2';
else
    disp('WARNING: NaN');
    c_hat = ckhat';
end

% Output
elapsed_time = toc();
[t, c] = sort_t_c(t_hat', c_hat);
sigma = se;

end

% compute fitness function f
function L = f(t, c, sigmah, n, y)
z = calc_z(t, c, sigmah, n);
L = -(y-z)*(y-z)';
end

% compute z
function z = calc_z(t, c, sigmah, n)
z = zeros(1, length(n));
for i = 1 : length(n)
    for k = 1 : length(t)
        z(i) = z(i) + c(k) * exp(-(n(i)-t(k))^2/(2*sigmah^2));
    end
end
end