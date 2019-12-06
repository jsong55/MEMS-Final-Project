function myfun
global gamma eta eta_y V_dc beta omega
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
eta_y = I_m*e*b*V_dc^2/(2*K_t*d^2*m*alpha_max);

myfun = @twodof;
x0 = [0.06,0.0125];
x = fsolve(myfun,x0)

end
function F = twodof(x)
global gamma eta eta_y beta omega
F(1) = x(1)-eta/(x(1)^2)*((1-x(2))/(1-x(2)-beta*x(1))-(1-x(2))/(1-x(2)-gamma*x(1))+log((1-x(2)-beta*x(1))/(1-x(2)-gamma*x(1))));
F(2) = (omega^2)*x(2)-(eta_y/x(1))*(1/(1-x(2)-beta*x(1))-1/(1-x(2)-gamma*x(1)));
end


