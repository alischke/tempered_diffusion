function y = evaluate_tempered_steady_state(x,alpha,lambda)
% y = evaluate_tempered_steady_state(x,alpha,lambda)

if (alpha > 2 || alpha <= 1)
    error('alpha must be in interval (1,2]');
end
if (lambda < 0)
    error('lambda must be non-negative')
end


%Map domain [L,R] to unit interval [0,1]
L = min(x);
R = max(x);
diam = R - L;
x1 = (x - L)/diam;
alpham1 = alpha -1;
alpham2 = alpha -2;

if (lambda > 0)

fac1 = mlf_star(alpha,alpham2,lambda.*x1);
fac2 = exp(-lambda.*x1);
u_steady1 = fac1 .* fac2;

fac1 = mlf_star(alpha,alpham1,lambda.*x1);
fac2 = exp(-lambda.*x1);
u_steady2 = fac1 .* fac2;

u_steady = -u_steady2 + u_steady1;

%Normalization constant
Aconst = exp(-lambda)*lambda^(alpha-2) * mlf(alpha,alpha,lambda^alpha);
y = u_steady./Aconst;

else
    y = alpham1.*x1.^(alpha-2);
end


y = y./diam;         %normalization constant

end



function y = mlf_star(alpha,beta,x)
% Evaluates modified Mittag-Leffler function E^{*}_{a, b} (x)
% See page 89 of Harish's Ph.D. thesis

arg = x.^alpha;
y = x.^beta .* mlf(alpha,beta + 1,arg);

end
