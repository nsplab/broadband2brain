function lnp = loglikelihood(ck, tk, sigmae, y, n)
y=y(:);
n=n(:);
ck=ck(:);
tk=tk(:);
N=length(y);
sigmah=1;
K = length(ck);
yprime=zeros(N,1);
for k = 1:K
   yprime=yprime+ck(k)*gausskernel(n-tk(k),sigmah);
end

lnp=-(N+1)*log(sigmae)-1/(2*sigmae^2)*sum((y-yprime).^2);

lnp=-lnp;