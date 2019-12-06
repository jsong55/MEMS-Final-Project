clear all
close all
K_t = 1.54*10^(-9);
a1 = 3*10^(-6);
a2 = 42*10^(-6);
a3 = 50*10^(-6);
d = 2.75*10^(-6);
b = 1000*10^(-6);
t = 1.5*10^(-6);
w = 1.55*10^(-6);
l = 65*10^(-6);
k_y = 6.49;
I_m = 2.5*10^(-20);
V_dc = 5.85; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 4.3*10^(-11);
alpha_max = d/a3;
gamma = a1/a3;
beta = a2/a3;
e = 8.85*10^(-12); %epsilon
eta = e*b*V_dc^2/(2*alpha_max^3*K_t);
omega_y = sqrt(k_y/m);
omega_t = sqrt(K_t/I_m);
omega = omega_y/omega_t;
x1 = 1;
eta_y = I_m*e*b*V_dc^2/(2*K_t*d^2*m*alpha_max);
%%
delta = 0.0125;
theta = linspace(-1,1.15);
V_1 = theta.^2/2+eta*log((1-beta*theta)./(1-gamma.*theta))./theta; %1dof
V_2 = theta.^2/2+eta*log((1-delta-beta*theta)./(1-delta-gamma.*theta))./theta; %2dof
%%
H1_1dof = -0.466;
H1_2dof = -0.466;
theta1_positive_1dof = sqrt(2*(H1_1dof-V_1));
theta1_negative_1dof = -sqrt(2*(H1_1dof-V_1));
theta1_positive_2dof= sqrt(2*(H1_2dof-V_2));
theta1_negative_2dof = -sqrt(2*(H1_2dof-V_2));
H2_1dof = -0.6;
H2_2dof = -0.6;
theta2_positive_1dof = sqrt(2*(H2_1dof-V_1));
theta2_negative_1dof = -sqrt(2*(H2_1dof-V_1));
theta2_positive_2dof= sqrt(2*(H2_2dof-V_2));
theta2_negative_2dof = -sqrt(2*(H2_2dof-V_2));
H3_1dof = -0.2;
H3_2dof = -0.2;
theta3_positive_1dof = sqrt(2*(H3_1dof-V_1));
theta3_negative_1dof = -sqrt(2*(H3_1dof-V_1));
theta3_positive_2dof= sqrt(2*(H3_2dof-V_2));
theta3_negative_2dof = -sqrt(2*(H3_2dof-V_2));
% H4_1dof = 0;
% H4_2dof = 0;
% theta4_positive_1dof = sqrt(2*(H4_1dof-V_1));
% theta4_negative_1dof = -sqrt(2*(H4_1dof-V_1));
% theta4_positive_2dof= sqrt(2*(H4_2dof-V_2));
% theta4_negative_2dof = -sqrt(2*(H4_2dof-V_2));
H5_2dof = -0.4822;
theta5_positive_2dof= sqrt(2*(H5_2dof-V_2));
theta5_negative_2dof = -sqrt(2*(H5_2dof-V_2));
%%
figure(2)
subplot(2,1,1)
plot(theta,V_1)
grid on
hold on
yline(0.-0.466,':.r');
yline(-0.6,'-.b');
yline(-0.2,'-.g');
xlabel('\theta')
ylabel('V(\theta)')
axis on
hold off
subplot(2,1,2)
plot(theta,theta1_positive_1dof,'r')
ylabel('$\dot{\theta}$','interpreter','latex')
xlabel('\theta')
grid on
hold on
plot(theta,theta1_negative_1dof,'r')
hold on
plot(theta,theta2_positive_1dof,'b')
hold on
plot(theta,theta2_negative_1dof,'b')
hold on
plot(theta,theta3_positive_1dof,'g')
hold on
plot(theta,theta3_negative_1dof,'g')
% hold on
% plot(theta,theta4_positive_1dof,'k')
% hold on
% plot(theta,theta4_negative_1dof,'k')
hold off
%%
figure(3)
subplot(2,1,1)
plot(theta,V_2)
grid on
hold on
yline(-0.466,'-.r');
yline(-0.4822,'-.b');
yline(-0.6,'-.g');
yline(-0.2,'-.k');
xlabel('\theta')
ylabel('V(\theta)')
axis on
hold off
subplot(2,1,2)
plot(theta,theta1_positive_2dof,':m')
grid on
hold on
plot(theta,theta1_negative_2dof,':m')
hold on
plot(theta,theta2_positive_2dof,'g')
hold on
plot(theta,theta2_negative_2dof,'g')
hold on
plot(theta,theta3_positive_2dof,'b')
hold on
plot(theta,theta3_negative_2dof,'b')
% hold on
% plot(theta,theta4_positive_2dof,'k')
% hold on
% plot(theta,theta4_negative_2dof,'k')
% hold on
plot(theta,theta5_positive_2dof,'r')
hold on
plot(theta,theta5_negative_2dof,'r')
ylabel('$\dot{\theta}$','interpreter','latex')
xlabel('\theta')
hold off
