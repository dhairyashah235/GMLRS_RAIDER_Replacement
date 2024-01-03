clc; clear
cd = 0.96;
g = 1.18; 
e = 6.2;
pc = (113.7); % lbf/in2
p0 = 0.020; % Pressure in psi
gamma0 = 1.4;
ae = 0.44; % nozzle exit area in in2
at = 0.077; % nozzle throat area in in2
a = 118*12; % in/s
rho = 0.000119; %density at cruise psi
syms m
x = (ae/at) == (((g+1)/2)^(-(g+1)/(2*(g-1))))*(((1+(((g-1)/2)*m^2))^((g+1)/(2*(g-1))))/m);
var = isolate(x,m);
Me = vpasolve(var,m);
% exit pressure calculations
q = 0.5*rho*(Me*a)^2;
ps = p0;
pt = (q+ps);
pe = pt*(1+((g-1)/2)*Me^2)^(-g/(g-1));

n1 = (2*g^2)/(g-1);
n2 = (2/(g+1))^((g+1)/(g-1));
n3 = 1-((pe/pc)^((g-1)/g));
T = cd*pc*at*((n1*n2*n3)^0.5 + (pe/pc)*e - (p0/pc)*e);
T=T*32.3

% % Specific Impulse calculations
mdot = 7.84; %lbm/s

G = (g+1)/(g-1);
R= 8.314;
T4 = 3232.78;
% c = 3.28084*sqrt((R*T4/g)*((g+1)/2)^G);
c = 5200; %ft/s
gr = 32.2;

Isp = cd*(((n1*n2*n3)^0.5)+(pe/pc)*e-(p0/pc)*e)*c/gr
