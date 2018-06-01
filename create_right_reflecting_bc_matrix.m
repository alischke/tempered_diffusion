function bmat = create_right_reflecting_bc_matrix(diam,n,lambda,alpha,model)
% Create the iteration matrix for a one-sided, negative tempered fractional 
% diffusion equation assuming reflecting BCs.
% n = number of gridpoints
% 1 < alpha <= 2 
% model is either 'norm' for 2-term normalized TFD eqn or 'cent' for 3-term
% centered normalized TFD eqn


if (alpha >2 || alpha <= 1)
    error('alpha must be greater than one and less than or equal to two')
end

% switch between normalized and centered normalized models
if (strcmp(model,'norm'))
    cflg = 0;
elseif (strcmp(model,'cent'))
    cflg = 1;
end

h = diam/n;
lh = lambda*h;
alpham1 = alpha-1;

np1 = n +1;
gwgts = create_grunwald_weights(alpha,n);

kappa = zeros(1,n+1);
for i = 0:n
    i1 = i+1;
    kappa(i1) = sum(gwgts(1:i1)'.*exp(-lh*(-1:i-1)));
end

bmat = zeros(np1);
for i = 0:n
    for j = 1:n
        i1 = i + 1;
        j1 = j + 1;
        k = i-j+1;
        k1 = k+1;
        
        if (i >= j-1)
            bmat(i1,j1) = gwgts(k1)*exp(-lh*(i-j));
        end
        
        if (i == j)
            bmat(i1,j1) = gwgts(k1) - exp(lh)*(1-exp(-lh))^alpha - ...
                cflg*alpha*(lh)^alpham1;
        end
        
        if (i == j+1)
            bmat(i1,j1) = gwgts(k1)*exp(-lh*(i-j)) + cflg*alpha*lh^alpham1;
        end
    end
    j = 0;
    i1 = i+1;
    j1 = j+1;
    bmat(i1,j1) = exp(lh)*(1-exp(-lh))^alpha - kappa(i1);
end

% special cases
i = 0; j = 0;
i1 = i+1; j1 = j+1;
bmat(i1,j1) = -gwgts(1)*exp(lh);

% i = 1; j = 0;
% i1 = i+1; j1 = j+1;
% bmat(i1,j1) = gwgts(3)*exp(-lh) + alpha*lh^alpham1;

i = n; j = n-1;
i1 = i+1; j1 = j+1;
bmat(i1,j1) = gwgts(3)*exp(-lh) + cflg*alpha*lh^alpham1;

i = n; j = n;
i1 = i+1; j1 = j+1;
bmat(i1,j1) = gwgts(2) - exp(lh)*(1-exp(-lh))^alpha - cflg*alpha*lh^alpham1 + ...
    gwgts(1)*exp(lh);

end