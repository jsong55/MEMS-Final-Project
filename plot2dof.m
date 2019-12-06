x1 = [0.0015,0.0061,0.0139,0.0252,0.0404,0.06,0.0851,0.1174,0.16,0.2208,0.33749,0.3892];
x3 = [0.0003,0.0013,0.0029,0.0053,0.0084,0.0125,0.0176,0.0241,0.0324,0.044,0.0649,0.0686];
V_2 = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6];
figure(1)
subplot(1,2,1)
plot(V_2,x1)
grid on
xlabel('V_{dc}')
ylabel('\theta')
subplot(1,2,2)
plot(V_2,x3)
xlabel('V_{dc}')
ylabel('\delta')
grid on
%%
eta_t = e*b*V_2.^2/(2*alpha_max^3*K_t);
torque1dof = eta_t./x1.^2.*(1./(1-beta.*x1)-1./(1-gamma.*x1)+log((1-beta.*x1)./(1-gamma.*x1)));
torque2dof = eta_t./x1.^2.*((1-x3)./(1-x3-beta.*x1)-(1-x3)./(1-x3-gamma.*x1)+log((1-x3-beta.*x1)./(1-x3-gamma.*x1)));
plot(V_2,torque1dof);
grid on
hold on
plot(V_2,torque2dof)
legend('One-DOF model torque','Two-DOF model torque')
xlabel('V_{dc}')
ylabel('Nm')
hold off
%%
% eta_w = e*b*V_2.^2/(2*alpha_max^3*K_t);
% w_non=real(sqrt(-1+eta_w./(x1.^2).*((1./(1-x1).^2-1./(1-x1)))-2.*eta_w./(x1.^3).*(-1+(1./(1-x1))+log(1-x1))));
% fn=w_non/(2*pi)*sqrt(K_t/I_m);
% plot(V_2,w_non)
%%
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
V_dc = 3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 4.3*10^(-11);
alpha_max = d/a3;
gamma = a1/a3;
beta = a2/a3;
e = 8.85*10^(-12); %epsilon
eta = e*b*V_dc^2/(2*alpha_max^3*K_t);
omega_y = sqrt(k_y/m);
omega_t = sqrt(K_t/I_m);
omega = omega_y/omega_t;
%x1 = 1;
eta_y = I_m*e*b*V_dc^2/(2*K_t*d^2*m*alpha_max);
%%
delta = 0.0125;
theta = linspace(-1,1.15);
V_1 = theta.^2/2+eta*log((1-beta*theta)./(1-gamma.*theta))./theta; %1dof
V_2 = theta.^2/2+eta*log((1-delta-beta*theta)./(1-delta-gamma.*theta))./theta; %2dof
%%
H1_1dof = 0.251;
H1_2dof = 0.251;
theta1_positive_1dof = sqrt(2*(H1_1dof-V_1));
theta1_negative_1dof = -sqrt(2*(H1_1dof-V_1));
theta1_positive_2dof= sqrt(2*(H1_2dof-V_2));
theta1_negative_2dof = -sqrt(2*(H1_2dof-V_2));
H2_1dof = 0.35;
H2_2dof = 0.35;
theta2_positive_1dof = sqrt(2*(H2_1dof-V_1));
theta2_negative_1dof = -sqrt(2*(H2_1dof-V_1));
theta2_positive_2dof= sqrt(2*(H2_2dof-V_2));
theta2_negative_2dof = -sqrt(2*(H2_2dof-V_2));
H3_1dof = 0.1;
H3_2dof = 0.1;
theta3_positive_1dof = sqrt(2*(H3_1dof-V_1));
theta3_negative_1dof = -sqrt(2*(H3_1dof-V_1));
theta3_positive_2dof= sqrt(2*(H3_2dof-V_2));
theta3_negative_2dof = -sqrt(2*(H3_2dof-V_2));
H4_1dof = 0;
H4_2dof = 0;
theta4_positive_1dof = sqrt(2*(H4_1dof-V_1));
theta4_negative_1dof = -sqrt(2*(H4_1dof-V_1));
theta4_positive_2dof= sqrt(2*(H4_2dof-V_2));
theta4_negative_2dof = -sqrt(2*(H4_2dof-V_2));
H5_2dof = 0.2318;
theta5_positive_2dof= sqrt(2*(H5_2dof-V_2));
theta5_negative_2dof = -sqrt(2*(H5_2dof-V_2));
%%
figure(2)
subplot(2,1,1)
plot(theta,V_1)
grid on
hold on
yline(0.2506,':.r');
yline(0.1,'-.b');
yline(0.35,'-.g');
yline(0,'-.k');
yline(-0.1244,'-.r');
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
plot(theta,theta2_positive_1dof,'g')
hold on
plot(theta,theta2_negative_1dof,'g')
hold on
plot(theta,theta3_positive_1dof,'b')
hold on
plot(theta,theta3_negative_1dof,'b')
hold on
plot(theta,theta4_positive_1dof,'k')
hold on
plot(theta,theta4_negative_1dof,'k')
hold off
%%
figure(3)
subplot(2,1,1)
plot(theta,V_2)
grid on
hold on
yline(0.2319,'-.r');
yline(0.1,'-.b');
yline(0.35,'-.g');
yline(0,'-.k');
yline(-0.1244,'-.r');
yline(0.2506,':.m');
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
hold on
plot(theta,theta4_positive_2dof,'k')
hold on
plot(theta,theta4_negative_2dof,'k')
hold on
plot(theta,theta5_positive_2dof,'r')
hold on
plot(theta,theta5_negative_2dof,'r')
ylabel('$\dot{\theta}$','interpreter','latex')
xlabel('\theta')
hold off
