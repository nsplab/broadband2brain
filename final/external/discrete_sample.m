function x=discrete_sample(N,values,proba);
% DISCRETE_SAMPLE Samples discrete values with given probabilities
%
% x=discrete_sample(N,values,proba)
%
% Samples from {x1,...,xn} with probabilities {p1,...,pn}
% N: number of samples
% values: {x1,...,xn} 
% proba: {p1,...,pn}
%
% Author: C. FEVOTTE 
% Email: cf269 AT cam DOT ac DOT uk 

x=zeros(1,N);

% Cumulative distribution
distrib=cumsum(proba);

for k=1:N
    u=rand(1);
    j=1; while u >distrib(j); j=j+1; end;
    x(k)=values(j);
end

