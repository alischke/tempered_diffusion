function p = tempered_stable(x,t,alpha,lambda,x0,model)
% Evaluates fundamental solution for tempered stable diffusion
% This code depends on the stbl MATLAB library by Mark Veillette
% Download the stbl library here: https://github.com/markveillette/stbl
% x0 = location of initial impulse
% model = 'norm' is the normalized tempered fractional diffusion equation
%       = 'cent' is the centered normalized tempered fractional diffusion equation     
% James F. Kelly
% 12 May 2018
% Modified by Anna Lischke
% 28 May 2018

skewfac = abs(cos(pi*alpha/2)).^(1/alpha);
scalefac = t^(-1/alpha);
alpham1 = alpha -1;


fac1 = exp(-lambda*(x-x0));
if (strcmp(model,'norm'))
  fac2 = exp(-t*lambda^alpha);
  z = (x-x0).*scalefac;
elseif (strcmp(model,'cent'))
  fac2 = exp(t*alpham1*lambda^alpha);
  z = (x-x0 - t*alpha*lambda^alpham1).*scalefac;
end  

fac3 = scalefac.*stblpdf(z,alpha,1,skewfac,0);
p = fac1 .* fac2 .* fac3;


end
