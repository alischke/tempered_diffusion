function gwgts = create_grunwald_weights(alpha,n)

np1 = n +1;
gwgts = zeros(np1,1);
gwgts(1) = 1;
gwgts(2) = -alpha;
alpham1 = alpha -1;
for i = 2:n
    i1 = i +1;
    fac1 = - (alpha - i + 1)/i;
    gwgts(i1) = gwgts(i) * fac1;
end

end