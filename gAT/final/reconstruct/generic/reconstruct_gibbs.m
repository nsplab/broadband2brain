function [t, c, sigma, elapsed_time, t_inter, c_inter, sig_inter gibbs_iterTime] = reconstruct_gibbs(y, T, args, prev_data)
%function [t_hat, c_hat, sigma_hat, time_y, y, c_hat_orig] = FRI_gibbs(time, x_real, N, K, numiter, mean_iter, sigmah)
% Adapted from improvedgibbsfunc.m from "Estimating Signals with Finite Rate of Innovation from Noisy Samples", V. Y. F. Tan
% Most hard-coded values have been replaced but se (sigma estimate) is still initialized to a hard-coded constant

% This code has been adapted to fit generic framework
% ** indicates changes

% Uses gaussian kernel with samples at endpoints of intervals (and equally
% spaced in between)

% Inputs:
% time - time vector for x_real ** time interval is now [0, T)
% x_real - real signal
% N - number of samples
% K - number of spikes
% numiter - number of iterations
% mean_iter - parameters estimated by taking mean of last mean_iter iterations, must have mean_iter <= numiter
% sigmah - width of gaussian kernel used to take samples

% Outputs:
% t_hat - spike times
% c_hat - spike amplitudes
% sigma_hat - estimate of sigma for gaussian noise
% c_hat_orig - original c_hat (before improved by least squares)

tic;

sigmah = args.sigmah;
N = args.L;
K = args.K;
numiter = args.numiter;
mean_iter = args.mean_iter;

t0 = 0; %min(time);  % Take min(time) to be 0
tf = T; %max(time) - min(time);
dt = get_dt(); %tf / (length(time) - 1);
%time = time - t0;
n = linspace(0, tf, N);  % Vector of nT values for n = 0, 1 ... N-1 (points in time at which to take samples y)

% Calculate y from x and h: y = x * h (convolution) ** samples are now taken separately and passed as argument
% h(t) = exp(-t^2 / 2*sigmah^2)
%y=zeros(length(n), 1);
%for j = 1 : length(n),
%    for k = 1 : length(x_real),
%        y(j) = y(j) + x_real(k) * exp(-(n(j)-time(k))^2 / (2 * sigmah^2)) * dt;
%    end
%end

%{
figure;
plot(y);
title('my z');
%}

% TESTING: add noise to y
%y = y + normrnd(0, 0.001, length(y), 1);

%{
figure;
plot(y);
title('samples y');
%}

% subplot(211); stem(n, z); title('(a) Clean signal z[n]');

% subplot(212); stem(n, y); title('(b) Noisy signal y[n]=z[n]+e[n]');



% c_0=0; c_1=0; c_2=0;

% t_0=0; t_1=0; t_2=0;

% Initial conditions

c=zeros(1,K);

% c=[6 12];

t=zeros(1,K);

% t=[10.5 13];

se=0.1;

% TESTING
%se = 0.001;

thetaest=[c t se];

% Discrete values of t used to draw from distribution
% 100 is number of values
t_sampled=linspace(0,tf,100);  % Edited by Alex Wein, changed 20 to tf



y=y';

% numiter=200;

%scalefactor=10^(-.5);

scalefactor=1;

% Save intermediate results
c_inter = [];
t_inter = [];
sig_inter = [];


for iter = 1:numiter
    
    iterID = tic;


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


        % Use discrete sampling

        %         for i = 1:length(t_sampled)

        %             lnprobs_t1(i)=-1/(2*se^2)*sum(gamma_0*exp(-((n-t_sampled(i)).^2)./(sigmah^2))+ ...

        %                 delta.*exp(-((n-t_sampled(i)).^2)./(2*sigmah^2)));

        %         end



        probs_t=exp(scalefactor*lnprobs_t)/sum(exp(scalefactor*lnprobs_t));  % probability of each t value in t_sampled

        %         figure; plot(t_sampled,probs_t)



        % Makes it faster

        %         if(iter>50)

        %             t_sampled_nonzero=t_sampled(probs_t>1e-16);

        %             probs_t_nonzeros=probs_t(probs_t>1e-16);



        % Implement rejection sampling here

        %             [m,maxindex]=max(probs_t_nonzeros);

        %             q=normpdf(t_sampled_nonzero, t_sampled_nonzero(maxindex), .3);

        %             q=q*m/max(q);

        %

        %             r=rand(1);

        %             rp=0; rq=1;

        %             if(length(find(q'-probs_t_nonzeros<-1e-2))~=0)

        %                 error('cq is not greater than p');

        %             end

        %

        %             while(r>rp/rq)

        %                 q_sample=normrnd(t_sampled_nonzero(maxindex), .3, 1);

        %                 rq=interp1(t_sampled_nonzero, q, q_sample,'linear','extrap');

        %                 rp=interp1(t_sampled_nonzero, probs_t_nonzeros, q_sample,'linear','extrap');

        %             end



        %             t(k)=q_sample;

        %             t(k)=discrete_sample(1, t_sampled_nonzero, probs_t_nonzeros);

        %         else

      if(norm(probs_t)==0)

        t(k)=rand(1)*tf;  % Edited by Alex Wein, changed 20 to tf
        disp('Warning: t_k chosen randomly');

      else

            t_sampled_nonzero=t_sampled;

            probs_t_nonzeros=probs_t;

            t(k)=discrete_sample(1, t_sampled_nonzero, probs_t_nonzeros);  % sample t_k from discrete values

      end

        %         end



    end

    %t_inter = [t_inter t'];  % save intermediate results




    % Update sigmae

    d=zeros(1,length(y));



    for kprime=1:K

        d=(c(kprime))*exp(-((n-t(kprime)).^2)./(2.*sigmah^2))+d;

    end



    psi=N/2;   % phi in paper (eq 27)

    phi=2/sum((y - d).^2);  % 2 / sum_over_n{(y_n - s_n)^2}, which is 1/lambda in paper (eq 28)

    se=(gamrnd(psi,phi))^(-.5);

    
    %sig_inter = [sig_inter se];  % save intermediate results

    gibbs_iterTime = toc(iterID);

    thetaest(iter,:) = [c  t  se];


%     drawnow;

%     subplot(321); plot(iter,c, 'rs', 'MarkerSize',2); hold on;   plot(iter, ck, 'bo', 'MarkerSize',2);

%     subplot(322); plot(iter,t, 'rs','MarkerSize',2); hold on;   plot(iter, tk, 'bo','MarkerSize',2);

%     subplot(3,2,3); plot(iter,se, 'rs','MarkerSize',2); hold on; plot(iter, sigmaek, 'bo','MarkerSize',2);

%     subplot(3,2,4); plot(iter,loglikelihood(c, t, se, y, n) , 'rs','MarkerSize',2); hold on;

    %     plot(iter, optobjfun, 'bo','MarkerSize',2);



    ll(iter)=loglikelihood(c,t,se,y,n);

end



% subplot(321); plot(thetaest(:,1:K), 'r-' ); hold on; plot(repmat(ck,numiter,1), 'b-');

% axis([0 numiter 0 max(ck)+2]);

% title('(a) c_k');

% subplot(322); plot(thetaest(:,K+1:2*K), 'r-'); hold on; plot(repmat(tk,numiter,1), 'b-');

% axis([0 numiter 0 max(tk)+2]); title('(b) t_k');

% subplot(3,2,3); plot(thetaest(:,end), 'r-'); hold on; plot(repmat(sigmaek,numiter,1), 'b-');

% title('(c) \sigma_e');

% axis([0 numiter 0 sigmaek+2]);

% subplot(3,2,4); plot(ll , 'r-'); %hold on; plot(repmat(optobjfun,numiter,1), 'b-');

% title('(d) -log(p)');

% axis([0 numiter min(ll )-10 max(ll )+10]);

if mean_iter == 1
    mtheta = thetaest(end,:);
else
    mtheta=mean(thetaest(end-mean_iter+1:end,:));
end

% Original outputs
ckhat=mtheta(1:K);
tkhat=mtheta(K+1:2*K);
sigmahat=mtheta(end);

% New outputs
c_hat_orig = mtheta(1:K);  % c_hat is improved below
t_hat = mtheta(K+1:2*K) + t0;
sigma_hat = mtheta(end);
time_y = n + t0;

%{

time = linspace(-10,30,100*N);

zhat=zeros(length(time),1);

ztrue=zeros(length(time),1);

for k = 1:K

    ztrue=ztrue+ck(k)*gausskernel(time-tk(k),sigmah)';

    zhat=zhat+ckhat(k)*gausskernel(time-tkhat(k),sigmah)';

end

%}

% subplot(325); plot(time, ztrue, 'b', time, zhat, 'r--'); axis([0 24 0 20]);

% % title('Plot of reconstructed z(t) and actual z(t)');

% legend('z(t)', 'z^{hat}(t)');



%error1=(norm(zhat-ztrue)/norm(ztrue))^2;

% title(strcat('(e) MMSE error=',num2str(error1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TESTING
%c_hat = c_hat_orig;
%return;

Hmatrix=zeros(N,K);



for k = 1:K

    for nn=1:N

        Hmatrix(nn,k)=gaussian(tkhat(k)-n(nn), args);

    end

end


ckhat2=Hmatrix\y';


%{

zhat2=zeros(length(time),1);

for k = 1:K

    zhat2=zhat2+ckhat2(k)*gausskernel(time-tkhat(k),sigmah)';

end

%}


% subplot(3,2,6); plot(time, ztrue, 'b', time, zhat2, 'r--'); axis([0 24 0 20]);

% legend('z(t)', 'z^{hat}(t)');



%error2=(norm(zhat2-ztrue)/norm(ztrue))^2;

% title(strcat('(f) LLSE error=',num2str(error2)));


% Output improved c_hat
c_hat = ckhat2';

% Output
[t, c] = sort_t_c(t_hat', c_hat);
sigma = sigma_hat;
elapsed_time = toc();
%gibbs_time = elapsed_time