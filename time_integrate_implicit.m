function [usnap,t] = time_integrate_implicit(u0,src,bt,deltat,nt,tout)
% [usnap,t] = time_integrate_implicit(u0,bt,deltat,nt,tout)
% Implicit Euler time-integration for fractional diffusion equation
% 12 February
nx = length(u0);
nsnap = length(tout);
usnap = zeros(nx,nsnap);

%start the clock
t = 0.0;
u = u0;

%Create LHS matrix
a_lhs = eye(nx) - bt;



for it = 1:nt
%    it
    s = src(:,it+1);
    rhs = u+s;
    up = mldivide(a_lhs,rhs);
    t = t + deltat;
    u = up;
   
   
   %store snapshots
   for isnap = 1:nsnap
     if (abs(t - tout(isnap))<1e-10)
       usnap(:,isnap) = u;
     end
   end
   
   
   
end


end
