% Compute solution to tempered space-fractional diffusion equation on a bounded
% domain, normalized and centered normalized tempered fractional derivatives
% Anna Lischke 
% June 1, 2018

clear all;
close all;

momsflg = 0;                  % turns source term for moms on (1) or off (0)

alpha = 1.5;                  %fractional order
lambda = 1;

n = 100;                      
nx = n + 1;                   %number of grid points
p = 0;                        %weight: p =1 is positive FD, p = 0 is negative FD, p = 1/2 is 
                              %fractional Laplacian
Cdiff = 1.0;                  %diffusion coefficient
deltat = 1e-3;                %time step


model = 'norm';                 %'norm'=normalized TFD and 'cent'=centered normalized TFD

xleft = 0;                      %set up spatial grid
xright = 1;
diam = xright -xleft;
np1 = n + 1;
h = diam/n;                     %grid spacing

bc_type = 'rr';                 %'a' = absorbing and 'r' = reflecting.  Specify both endpoints
integrate_type = 'i';           %'i' = implicit Euler and 'e' = explicit Euler 
x = xleft + h.*(0:n)';          %grid

tout = [0.1 0.3 0.5 0.8 1];          %specify output times
nsnap = length(tout);
nt = ceil(tout(nsnap)./deltat);      %number of timesteps  
t = 0:deltat:nt*deltat;

% Initial Conditions
if momsflg == 1
    bta = 2;
    u0 = moms_ic(x,alpha,bta,lambda);
    src = source(x,t,alpha,bta,lambda); % define source function for moms
elseif momsflg == 0
	% tent initial condition
    u0 = tent(x);
    u0 = u0/(sum(u0).*h);
    src = zeros(length(x),length(t));
	% impulse initial condition
%     u0 =zeros(size(x));
%     u0(n/2 + 1) = 1/h;
end

ini_mass = sum(u0)*h;           %initial mass

bt = create_itmatrix_tempered(diam,p,Cdiff,deltat,h,n,alpha,lambda,bc_type,model);

% Time-stepping
if (strcmp(integrate_type,'e'))
  cfl = h^alpha / (Cdiff*alpha);

  if (deltat > cfl)
     error('time step is violating CFL limit')
     deltat
     cfl
  end
  [usnap, t] = time_integrate(u0,src,bt,deltat,nt,tout);
elseif (strcmp(integrate_type,'i'))  
  [usnap,t] = time_integrate_implicit(u0,src,bt,deltat,nt,tout);
else
   error('integrate_type must be i or e')
end

u = usnap(:,nsnap);
final_mass = sum(u)*h

%steady state solution for reflecting BCs, positive deriv (p = 1), model = 'norm'
if (p == 1) && strcmp(bc_type,'rr') && strcmp(model,'norm')
    u_steady = evaluate_tempered_steady_state(x,alpha,lambda);
else u_steady = NaN*ones(size(x));
end

%scale by mass
u_steady = ini_mass.*u_steady;

%Decimate for plotting

xd = x(1:8:nx);
us = u_steady(1:8:nx);
moms = exp(-tout(1))*u0(1:8:nx);

%Plot solutions

figure(1)
h1 = plot(x,u0,'-',x,usnap(:,1),':',x,usnap(:,2),'-.',x,usnap(:,3),'--',...
    x,usnap(:,4),'-',x,usnap(:,5),':');

set(h1,'LineWidth',3)
xlabel('x')
ylabel('u(x,t)')
title(['\alpha = ',num2str(alpha),',  \lambda = ',num2str(lambda),',  ',model])
grid on
set(gca,'FontSize',20)

%Plot solutions

figure(2)
h1 = plot(x,u0,'-',x,usnap(:,1),xd,moms,'o');

set(h1,'LineWidth',3)
set(gca,'FontSize',15)
xlabel('x')
ylabel('u(x,t)')
leg=legend('t = 0',['t = ',num2str(tout(1))],'MOMS');
set(leg,'Location','best')
title(['\alpha = ',num2str(alpha),',  \lambda = ',num2str(lambda),',  ',model],'FontSize',18)
grid on
set(gca,'xtick',[xleft:((xright-xleft)/10):xright])
set(gca,'FontSize',20)


% % Plot comparisons to tempered stable distribution (only for left-sided
% case)
% if p==1
%     x0 = xleft+diam/2;
%     for i = 1:length(tout)
%         tmpstbl(:,i) = tempered_stable(x,tout(i),alpha,lambda,x0,model);
%     end
% 
%     if strcmp(integrate_type,'i')
%         scheme = 'imp';
%     else scheme = 'exp';
%     end
% 
%     figure(3)
%     plot(x,tmpstbl(:,2),x,tmpstbl(:,3),x,tmpstbl(:,4),...
%         x,usnap(:,2),'--',x,usnap(:,3),'--',...
%         x,usnap(:,4),'--','LineWidth',2)
%     legend(['tmp stbl t = ',num2str(tout(2))],...
%         ['tmp stbl t = ',num2str(tout(3))],...
%         ['tmp stbl t = ',num2str(tout(4))],...
%         [scheme,' Euler t = ',num2str(tout(2))],...
%         [scheme,' Euler t = ',num2str(tout(3))],...
%         [scheme,' Euler t = ',num2str(tout(4))],'Location','best');
%     xlabel('x'); ylabel('u(x,t)');
%     title(['\alpha = ',num2str(alpha),',  \lambda = ',num2str(lambda),',  ',model],'FontSize',18)
%     grid on
%     set(gca,'FontSize',20)
% end








