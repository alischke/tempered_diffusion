function bt = create_itmatrix_tempered(p,Cdiff,deltat,h,n,alpha,lambda,bc_type,model)

% bc_type ='rr', 'ar', 'ra', 'aa'

nx = n +1;
q = 1 - p;
betap = p*Cdiff *h.^(-alpha) * deltat;
betaq = q*Cdiff *h.^(-alpha) * deltat;


%STEP 1: Create iteration matrix w/ reflecting BCs for a postive (left) FD

bmat_left = create_left_reflecting_bc_matrix(n,lambda,alpha,model);
% bmat_left = create_int_reflecting_bc_matrix(n,lambda,alpha);

%STEP 2: Use symmetry to create iteration matrix w/ reflecting BCs for a negative (right) FD
bmat_right = flipud(fliplr(bmat_left));

%Step 3: Zero columns if there are absorbing BCs
if (strcmp(bc_type(1),'a'))
    bmat_left(:,1) = 0;
    bmat_right(:,1) = 0;
end
if (strcmp(bc_type(2),'a'))
    bmat_left(:,nx) = 0;
    bmat_right(:,nx) = 0;
end

% Assign weights to each matrix
btp =betap.*transpose(bmat_left) ;
btq = betaq.*transpose(bmat_right);
bt = btp + btq;

end