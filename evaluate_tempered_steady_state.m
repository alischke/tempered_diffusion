function y = evaluate_tempered_steady_state(x,alpha,lambda)
% y = evaluate_tempered_steady_state(x,alpha,lambda)

if (alpha > 2 || alpha <= 1)
    error('alpha must be in interval (1,2]');
end
if (lambda < 0)
    error('lambda must be non-negative')
end
L = min(x);
R = max(x);
diam = R-L;
alpham1 = alpha -1;
alpham2 = alpha -2;

if (lambda > 0)

fac1 = mlf_star(alpha,alpham2,lambda.*x);
fac2 = exp(-lambda.*x);
u_steady1 = fac1 .* fac2;

fac1 = mlf_star(alpha,alpham1,lambda.*x);
fac2 = exp(-lambda.*x);
u_steady2 = fac1 .* fac2;

u_steady = -u_steady2 + u_steady1;

%Normalization constant

Aconst = (lambda^alpham2)*exp(-lambda*R)*(R^alpham1)*mlf(alpha,alpha,(lambda*R)^alpha);
y = u_steady./Aconst;

else
    y = alpham1.*x.^(alpha-2);
    y = y./diam;
end


end



function y = mlf_star(alpha,beta,x)
% Evaluates modified Mittag-Leffler function E^{*}_{a, b} (x)
% See page 89 of Harish's Ph.D. thesis

arg = x.^alpha;
y = x.^beta .* mlf(alpha,beta + 1,arg);

end
