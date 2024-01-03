clc;clear
Wto = 666; %lb
Wp = 201; % lb
a = 1125; % ft/s
V = a*4; % cruise speed
%Cl = 1.15; % lift coeff
Sref = 62.8; %ft2
%rho_cruise = (2*(Wto-0.5*Wp)/(V^2*Cl*Sref));
%h = 150000; % cruise altitude
P0 = 0.020; % Pressure in psi
gamma0 = 1.4;
T0 = 439.89; % Free stream temp (R)
Hf = 17900;% RJ-5 fuel heating value
M0 = 4; % free stream mach number
cp = 0.302; % specific heat of air at constant pressure in the combustor
fbya = 0.067; % fuel to air ratio
T4 = T0*(1+((gamma0-1)/2)*M0^2) + (Hf/cp)*(fbya);
temp_ratio = T4/T0;

Th = 235*32; % thrust in lbs

var = (temp_ratio/(1+((gamma0-1)/2)*M0^2));
A0inv = P0*gamma0*M0^2*(sqrt(var)-1)/Th;
A0 = 1/A0inv
rho = 0.000119; %density at cruise psf
mdota = rho*Sref*V/144

Ait = Sref*1.728*1.5*(1+(0.2*1.5^2))^-3

%% Section 4.4

M3tc = 0.461*((1+0.2*M0^2)/temp_ratio)^0.5
A3 = ((1.728*M3tc)/((1+0.2*M3tc^2)^3)*Ait)^-1
T3 = (1+0.2*M0^2)*T0;
a1 = 1.822*(T4/T3)*M3tc^2-1.175;
b1 = 2.7*(T4/T3)*M3tc^2;
c1 = (T4/T3)*M3tc^2;

M4 = ((-b1-(b1^2-(4*a1*c1))^0.5)/(2*a1))^0.5;
M4 = abs(M4);

% % Assuming gamma = 1.35
g2 = 1.35;
exp = g2/(g2-1);
pratio = ((1+((g2-1)/2)*M4^2)^(exp))/(1+g2*M4^2)

