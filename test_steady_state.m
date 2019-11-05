% Unit test for  evaluate_tempered_steady_state
% James F. Kelly
% 4 November 2019

clear all;
nx = 1000001;

alpha = rand + 1
lambda = 3*rand


L = 0; R = 1;
x = linspace(L,R,nx);
dx = x(2) - x(1);
y = evaluate_tempered_steady_state(x,alpha,lambda);
int1 = sum(y(2:nx)).*dx;

midp = 2.*(rand -1);
rad = 2.*rand;

L = midp - rad;
R = midp + rad;
x = linspace(L,R,nx);
dx = x(2) - x(1);
y = evaluate_tempered_steady_state(x,alpha,lambda);
int2 = sum(y(2:nx)).*dx;

figure(1)
plot(x,y)
axis([L R 0 10])
set(gca,'FontSize',20)
xlabel('x')
ylabel('steady state u_{\infty} (x)')
titlestr = ['\alpha = ' num2str(alpha) ' lambda = ' num2str(lambda)];
title(titlestr)



if max(abs(int1 - int2)) > 1e-10
    error('check scaling in evaluate_tempered_steady_state')
end
