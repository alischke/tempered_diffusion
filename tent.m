function y = tent(x)
% Modified May 25, 2018 for arbitrary interval [L,R].


nx = length(x);
L = min(x);
R = max(x);
diam = R - L;
x1 = (x - L)/diam;


y = zeros(size(x));
for ix = 1:nx
    xx = x1(ix);
    if (xx>0.3 && xx <= 0.5)
        y(ix) = 25*xx - 7.5;
    elseif (xx>0.5 && xx < 0.7)
        y(ix) = -25*xx + 17.5;
    end
end


y = y./diam;

end
