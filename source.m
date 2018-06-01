function s = source(x,t,alpha,bta,lambda)

u0 = moms_ic(x,alpha,bta,lambda);

fac1 = exp(-t)*(lambda^alpha-1);
fac2 = exp(-t-lambda*x);
fac3 = (((2^bta)/(1+bta)) - (1+x).^bta);

s = fac1.*u0 - fac2.*fac3;


end