function bmat = create_right_reflecting_bc_matrix(n,lambda,alpha)
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

h = 1/n;
lh = lambda*h;

np1 = n +1;
gwgts = create_grunwald_weights(alpha,n);
gwgtsm1 = create_grunwald_weights(alpha-1,n);

kappa = zeros(1,n+1);
for i = 0:n
    i1 = i+1;
    kappa(i1) = sum(gwgts(1:i1)'.*exp(-lh*(-1:i-1)));
end

bmat = zeros(np1);
for i = 0:n
    for j = 1:(n-1)
        i1 = i + 1;
        j1 = j + 1;
        
        if (i >= (j-1))
            k = i - j + 1;
            k1 = k + 1;
            bmat(i1,j1) = gwgts(k1)*exp(-lh*(i-j));
        end
        
        if i == j
            bmat(i1,j1) = gwgts(2) - exp(lh)*(1-exp(-lh))^alpha + ...
                cflg*alpha*lh^(alpha-1);
        end
        
        if i == j+1
            bmat(i1,j1) = gwgts(3)*exp(-lh) - cflg*alpha*lh^(alpha-1);
        end
            
    end
    if (i <= n)
    j = 0;
    i1 = i + 1;
    j1 = j + 1;
    bmat(i1,j1) = exp(lh)*(((1-exp(-lh))^alpha))-kappa(i1);
    end
end

% special cases - RL flux
i = 0; j = 0;
i1 = i+1; j1 = j+1;
bmat(i1,j1) = -gwgts(1)*exp(lh);

i = n; j = 0;
i1 = i+1; j1 = j+1;
k = i;
k1 = k+1;
bmat(i1,j1) = exp(lh)*(((1-exp(-lh))^alpha)) - kappa(k1);

i = n-1; j = n;
i1 = i+1; j1 = j+1;
bmat(i1,j1) = gwgts(1)*exp(lh);

i = n; j = n;
i1 = i+1; j1 = j+1;
k = j-i;
k1 = k+1;
np1 = n+1;
bmat(i1,j1) = gwgts(1)*exp(lh) - exp(lh)*((1-exp(-lh))^alpha)+gwgts(2)+...
    cflg*alpha*lh^(alpha-1);



end