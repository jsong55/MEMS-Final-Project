clear all
close all
syms x theta
K_t = 1.4308*10^(-5);
a1 = 430*10^(-6);
a2 = 580*10^(-6);
a3 = 600*10^(-6);
d = 3.42*10^(-6);
b = 1300*10^(-6);

alpha_max = d/a3;
gamma = a1/a3;
beta = a2/a3;
e = 8.85*10^(-12); %epsilon

x = linspace(0,0.3);
V_dc = sqrt((2*K_t*alpha_max^3*x.^3)./(e*b*((1./(1-beta*x))-(1./(1-gamma*x))+log((1-beta*x)./(1-gamma*x)))));
figure(1)
plot(x,V_dc)
grid on
xlabel('{\theta}')
ylabel('{V_{dc}}')
x1 = linspace(0.01,0.3961);
V_dc_1 = sqrt((2*K_t*alpha_max^3*x1.^3)./(e*b*((1./(1-beta*x1))-(1./(1-gamma*x1))+log((1-beta*x1)./(1-gamma*x1)))));
figure(2)
plot(V_dc_1,x1)
grid on
ylabel('{\theta}')
xlabel('{V_{dc}}')
x2 = linspace(0.3961,1);
V_dc_2 = sqrt((2*K_t*alpha_max^3*x2.^3)./(e*b*((1./(1-beta*x2))-(1./(1-gamma*x2))+log((1-beta*x2)./(1-gamma*x2)))));
figure(3)
plot(V_dc_2,x2,'--')
grid on
ylabel('{\theta}')
xlabel('{V_{dc}}')
%%
eta_w = x.^3./(1./(1-beta*x)-1./(1-gamma*x)+log((1-x)/(1-gamma*x)));
w_non=real(sqrt((-1+eta_w./(x.^2).*((1./(1-x).^2-1./(1-x))))-(2.*eta_w./(x.^3).*(-1+(1./(1-x))+log(1-x)))));
I_m=4.372*10^(-15);
fn=w_non/(2*pi)*sqrt(K_t/I_m);
plot(V_dc,fn)
%%
VDC=18;
eta = e*b*VDC^2/(2*alpha_max^3*K_t);
theta = linspace(-1,1);
V = theta.^2/2+eta*log((1-beta*theta)./(1-gamma.*theta))./theta;
figure(4)
subplot(2,1,1)
plot(theta,V)
grid on
hold on
yline(-0.1777,'-.b');
yline(-0.4,'-.b');
yline(0.2,'-.b');
yline(0,'-.b');
xlabel('\theta')
ylabel('V(\theta)')
axis on
% title('V_{dc}=10V')
hold off
%%
H1 = -0.1777;
theta1_dot_positive = sqrt(2*(H1-V));
theta1_dot_negative = -sqrt(2*(H1-V));
subplot(2,1,2)
plot(theta,theta1_dot_positive,'b')
axis on
grid on
hold on
plot(theta,theta1_dot_negative,'b')
hold on
H2 = -0.4;
theta2_dot_positive = sqrt(2*(H2-V));
plot(theta,theta2_dot_positive,'b')
hold on
theta2_dot_negative = -sqrt(2*(H2-V));
plot(theta,theta2_dot_negative,'b')
hold on
H3 = 0.2;
theta3_dot_positive = sqrt(2*(H3-V));
plot(theta,theta3_dot_positive,'b')
hold on
theta3_dot_negative = -sqrt(2*(H3-V));
plot(theta,theta3_dot_negative,'b')
hold on
H4 = 0;
theta4_dot_positive = sqrt(2*(H4-V));
plot(theta,theta4_dot_positive,'b')
hold on
theta4_dot_negative = -sqrt(2*(H4-V));
plot(theta,theta4_dot_negative,'b')
ylabel('$\dot{\theta}$','interpreter','latex')
xlabel('\theta')
title('V_{dc}=10V')
hold off
