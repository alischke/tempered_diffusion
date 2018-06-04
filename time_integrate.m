function [usnap,t] = time_integrate(u0,src,bt,deltat,nt,tout)

nx = length(u0);
nsnap = length(tout);
usnap = zeros(nx,nsnap);


%start the clock
t = 0.0;
u = u0;

for it = 1:nt
%    it
   s = src(:,it);
   rhs = bt * u;
   up = u + rhs + s*deltat;
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
