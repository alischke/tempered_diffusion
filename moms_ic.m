function u = moms_ic(x,alpha,bta,lambda)

p = ((1+x).^alpha)./gamma(1+alpha);

fac1 = exp(-lambda*x);
fac2 = (2^bta) / (1+bta);
fac3 = gamma(1+bta)/gamma(alpha+bta+1);

u = fac1.*(fac2.*p - fac3.*(1+x).^(alpha+bta));

end