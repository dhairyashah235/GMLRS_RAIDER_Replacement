%% Section 5
clc
close all
clear all

% Code Notes
% Code to solve for the AC loaction as a ratio of the distance to the body
% aerodynamic cenetr to the length of the nose of the missile.

% Variables
alpha = 0:1:90; % Angle of Attack Range (deg)
lB = 155; % Total Length of Body (in)
lN = 18.7; % Length of Nose (in)
RatioGMLRS = lB/lN % lB/lN of GMLRS (~)
Ratio1 = 1; % lB/lN of 1 (~)
Ratio2 = 2; % lB/lN of 2 (~)
Ratio5 = 5; % lB/lN of 5 (~)
Ratio10 = 10; % lB/lN of 10 (~)

% Equation
BACGMLRS = (0.63*(1-(sind(alpha)).^2))+(0.5*(RatioGMLRS)*(sind(alpha)).^2);
BAC1 = (0.63*(1-(sind(alpha)).^2))+(0.5*(Ratio1)*(sind(alpha)).^2);
BAC2 = (0.63*(1-(sind(alpha)).^2))+(0.5*(Ratio2)*(sind(alpha)).^2);
BAC5 = (0.63*(1-(sind(alpha)).^2))+(0.5*(Ratio5)*(sind(alpha)).^2);
BAC10 = (0.63*(1-(sind(alpha)).^2))+(0.5*(Ratio10)*(sind(alpha)).^2);

% Graph
figure
plot(alpha,BACGMLRS,'b')
hold on
plot(alpha,BAC1,'k')
plot(alpha,BAC2,'m')
plot(alpha,BAC5,'g')
plot(alpha,BAC10,'r')
hold off
ylim([0 5])
xlim([0 100])
yticks([0 0.63 1 2 3 4 5]);
xticks([0 10 20 30 40 50 60 70 80 90 100]);
ylabel([{'Distance to Body Aerodynamic Center/Length of Nose,' ...
    '{(x_{AC})_B}/l_N (~)'}]);
xlabel('Angle of Attack, \alpha (\circ)');
legend('l_B/l_N = 8.29 [GMLRS]','l_B/l_N = 1','l_B/l_N = 2', ...
    'l_B/l_N = 5','l_B/l_N = 10','Location','northwest');

%% Section 6
clc
clear
C_N_AoA_B = 2;%body'sNormal force coefficient dderivative contributors to pitching moment(/rad) 
X_CG = 77.5%Location of the Center of Gravity (inches)
L_N = 18.7;%Nose length (inches)
L_B = 157.15;%length of body-flare GMLRS rocket (inches)
d_F = 11.73;%flare diameter(inches)
d = 8.94;%missile diameter(inches)
X_F = 155;%loaction of the flare(inches)
L_F = 2.15;%flare length (inches)
AoA = [0:0.5:90];
%...................................Calculation For GMLRS Rocket........................
X_ac_b_over_L_N = 0.63*(1-(sind(AoA).^2))+(L_B/L_N)*(sind(AoA).^2)/2
X_ac_b = X_ac_b_over_L_N*L_N
C_N_AoA_F = 2*((d_F/d)^2-1);
X_ac_F = X_F+0.33*L_F*((d_F/d)*2+1)/((d_F/d)+1);
S_M = (C_N_AoA_F.*((X_CG-X_ac_F)/d)-C_N_AoA_B*(X_CG-X_ac_b)/d)/(C_N_AoA_B+C_N_AoA_F);
figure(1)
plot(AoA,S_M)
xlabel('Angle of Attack, \alpha(deg)','FontSize',26)
ylabel('Static Margin, S.M.(~)','FontSize',26)
set(gca,'FontSize',20)
set(gcf,'color','w')
X_ac = S_M*d+X_CG;
LB_LN_1  = 1;
LB_LN_2  = 2;
LB_LN_5  = 5;
LB_LN_10  = 10;
X_ac_b_over_L_N_1 = 0.63*(1-(sind(AoA).^2))+(LB_LN_1)*(sind(AoA).^2)/2;
X_ac_b_over_L_N_2 = 0.63*(1-(sind(AoA).^2))+(LB_LN_2)*(sind(AoA).^2)/2;
X_ac_b_over_L_N_5 = 0.63*(1-(sind(AoA).^2))+(LB_LN_5)*(sind(AoA).^2)/2;
X_ac_b_over_L_N_10 = 0.63*(1-(sind(AoA).^2))+(LB_LN_10)*(sind(AoA).^2)/2
figure(2)
K1 = plot(AoA,X_ac/L_N,'r-.')
hold on
K2 = plot(AoA,X_ac_b_over_L_N_1,'k-');
K3 = plot(AoA,X_ac_b_over_L_N_2,'g-');
K4 = plot(AoA,X_ac_b_over_L_N_5,'b-');
K5 = plot(AoA,X_ac_b_over_L_N_10,'r-');
K6 = legend([K1 K2 K3 K4 K5],{'GMLRS Rocket with Flare Stabilizer, total Length of Body/Length of Nose = 8.3','Total Length of Body/Length of Nose = 1','Total Length of Body/Length of Nose = 2', ...
    'Total Length of Body/Length of Nose = 5','Total Length of Body/Length of Nose = 10'},'Location','northwest');
hold off
set(gca,'FontSize',20)
set(gcf,'color','w')
ylim([0 5.2])
xlim([0 100])
grid on
xlabel('Angle of Attack, \alpha(deg)','FontSize',26);
ylabel('Distance to Body Aerodynamic center/Length of Nose, XAC/lN(~)','FontSize',20);

%% Section 7
% % Figure 2.25

clear all; 
clc;
%variables
Ssurf_t = 833.97/144; % in^2 surface platform area for tail 
Ssurf_w = 44.23/144; % in^2 surface platform area for wing/canard
M = 1:5; % Mach number
d = 8.94/12; % diameter in inches
Sref = pi/4*d^2; % Reference area in in^2
AR_w = 4.25;
AR_t = 12.8;
S_w = Sref/Ssurf_w; % reference area by wing surface area
S_t = Sref/Ssurf_t; %reference area by tail surface area

dCNbyda_sw = zeros (1,5); % dCNbyda for wing using slender wing theory
dCNbyda_lw = zeros (1,5); % dCNbyda for wing using linear wing theory
dCNbyda_st = zeros (1,5); % dCNbyda for tail using slender wing theory
dCNbyda_lt = zeros (1,5); % dCNbyda for tail using linear wing theory

for i = 1:5
    if  M(i) > sqrt(1+(8/(pi*AR_w))^2)   
        dCNbyda_lw(i) = (4/sqrt((M(i)^2)-1))*(Ssurf_w/Sref);
    elseif M(i) < sqrt(1+(8/(pi*AR_w))^2)
        dCNbyda_sw(i) = (pi*AR_w/2)*(Ssurf_w/Sref);
    end
end    

for i =1:5
    if M(i) > sqrt(1+(8/(pi*AR_t))^2)
        dCNbyda_lt(i) = (4/sqrt((M(i)^2)-1))*(Ssurf_t/Sref);
    elseif M(i) < sqrt(1+(8/(pi*AR_t))^2)
        dCNbyda_st(i) = (pi*AR_t/2)*(Ssurf_t/Sref);
    end
end
dCNbyda_lw = dCNbyda_lw*(Sref/Ssurf_w);
dCNbyda_sw = dCNbyda_sw*(Sref/Ssurf_w);
dCNbyda_lt = dCNbyda_lt*(Sref/Ssurf_t);
dCNbyda_st = dCNbyda_st(:,1)*(Sref/Ssurf_t);

    figure 
    grid on
    plot(M,dCNbyda_lw, "--black");
    hold on
    yline(dCNbyda_sw(:,1),"red","AR = 4.25")
    hold on
    plot(M,dCNbyda_lt, "--black");
    hold on
    yline(dCNbyda_st(:,1),"red","AR = 12.8")
    hold off
    xlim([0.1 5])
    
 xlabel('M, Mach Number')
 ylabel(' Non-dimensional Normal Force Coefficient Slope with Angle of Attack, per rad')
 legend('Linear Wing Thoery','Slender Wing Theory') 
   
% % Figure 2.26 
clear all;
clc;
% Variables
% Mach number range
M1 = 1.35;
M2 = 2;
M3 = 5;
Mach = [M1 M2 M3];
CN_w = zeros(3,91); % Coefficient of normal force for wing 
CN_t = zeros(3,91);% Coefficient of normal force for tail 
AOA = 0:90; % alpha dash is the same as AOA
AR_w = 4.25; % aspect ratio of the wing
AR_t = 12.79; % aspect ratio of the tail
Ssurf_w = 833.97; % wing surface platform area in in^2
Ssurf_t = 44.23; % tail curface planform area in in^2
d = 8.94; % diameter in inches
Sref = pi/4*d^2; % Reference area in in^2
S_w = Sref/Ssurf_w; % reference area by wing surface area
S_t = Sref/Ssurf_t; %reference area by tail surface area

for j = 1:91
    for i = 1:3
        if Mach(i) <= 1.35
            CN_w(i,j) = (((pi*AR_w/2)*abs(sind(AOA(j))*cosd(AOA(j))))+2*((sind(AOA(j)))^2));
            CN_t(i,j) = (((pi*AR_t/2)*abs(sind(AOA(j))*cosd(AOA(j))))+2*((sind(AOA(j)))^2));

        elseif Mach(i) == 2
            CN_w(i,j) = (((4*abs(sind(AOA(j))*cosd(AOA(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(AOA(j)))^2));
            CN_t(i,j) = (((4*abs(sind(AOA(j))*cosd(AOA(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(AOA(j)))^2));
        elseif Mach(i) == 5
            CN_w(i,j) = (((4*abs(sind(AOA(j))*cosd(AOA(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(AOA(j)))^2));
            CN_t(i,j) = (((4*abs(sind(AOA(j))*cosd(AOA(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(AOA(j)))^2));
        end
    end

end

figure 
plot(AOA,CN_w(1,:),"-.magenta")
hold on
plot (AOA,CN_w(2,:), "green")
hold on
plot (AOA, CN_w(3,:), "--blue")
grid on
hold off 
ylim([0 6]);
xlabel('Wing Effective Angle of Attack, (\alpha), Deg')
ylabel({'Wing Normal Force Coefficient for rocket baseline missile';'(C_N)_{wing} S_{ref} / S_w (~)'})
legend('M < 1.35, based on Slender Wing Theory + Newtonian Impact Theory','M = 2, based on Linear Wing Theory + Newtonion Impact Theory','M = 5, based on Linear Wing Theory + Newtonion Impact Theory') 

figure 
plot(AOA,CN_t(1,:),"-.magenta")
hold on
plot (AOA,CN_t(2,:), "green")
hold on
plot (AOA, CN_t(3,:),"--blue")
grid on
hold off
ylim([0 20]);
xlabel('Tail Effective Angle of Attack, (\alpha), Deg')
ylabel({'Tail Normal Force Coefficient for rocket baseline missile';'(C_N)_{tail} S_{ref} / S_t (~)'})
legend('M < 1.35, based on Slender Wing Theory + Newtonian Impact Theory','M = 2, based on Linear Wing Theory + Newtonion Impact Theory','M = 5, based on Linear Wing Theory + Newtonion Impact Theory') 

% Regenerating figure 2.26 for given deflection 
del_w = 5; % deflection in deg  
del_t = -2; % deflection in deg for tail
alpha_w = AOA+del_w; % effective angle of attack for wing/canard
alpha_t = AOA+del_t; % effective angle of attack for tail
for j = 1:91
    for i = 1:3
        if Mach(i) <= 1.35
            CN_w(i,j) = (((pi*AR_w/2)*abs(sind(alpha_w(j))*cosd(alpha_w(j))))+2*((sind(alpha_w(j)))^2));
            CN_t(i,j) = (((pi*AR_t/2)*abs(sind(alpha_t(j))*cosd(alpha_t(j))))+2*((sind(alpha_t(j)))^2));

        elseif Mach(i) == 2
            CN_w(i,j) = (((4*abs(sind(alpha_w(j))*cosd(alpha_w(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(alpha_w(j)))^2));
            CN_t(i,j) = (((4*abs(sind(alpha_t(j))*cosd(alpha_t(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(alpha_t(j)))^2));
        elseif Mach(i) == 5
            CN_w(i,j) = (((4*abs(sind(alpha_w(j))*cosd(alpha_w(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(alpha_w(j)))^2));
            CN_t(i,j) = (((4*abs(sind(alpha_t(j))*cosd(alpha_t(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(alpha_t(j)))^2));
        end
    end

end

figure 
plot(alpha_w,CN_w(1,:),'-.magenta')
hold on
plot (alpha_w,CN_w(2,:),'green')
hold on
plot (alpha_w, CN_w(3,:),'--blue')
grid on
hold off 
xlim([0 90]);
ylim([0 6]);
hold off
xlabel('Wing Effective Angle of Attack, (\alpha), Deg')
ylabel({'Wing Normal Force Coefficient for GMLRS';'(C_N)_{wing} S_{ref} / S_w (~)'})
legend('M < 1.35, based on Slender Wing Theory + Newtonian Impact Theory','M = 2, based on Linear Wing Theory + Newtonion Impact Theory','M = 5, based on Linear Wing Theory + Newtonion Impact Theory') 

figure 
plot(alpha_t,CN_t(1,:),"-.magenta")
hold on
plot (alpha_t,CN_t(2,:),"green")
hold on
plot (alpha_t, CN_t(3,:),"--blue")
grid on
hold off
xlim([0 90]);
ylim([0 20]);

xlabel('Tail Effective Angle of Attack, (\alpha), Deg')
ylabel({'Wing Normal Force Coefficient for GMLRS';'(C_N)_{tail} S_{ref} / S_t (~)'})
legend('M < 1.35, based on Slender Wing Theory + Newtonian Impact Theory','M = 2, based on Linear Wing Theory + Newtonion Impact Theory','M = 5, based on Linear Wing Theory + Newtonion Impact Theory')

%% Section 8
clc
clear all

%Define Aspect Ratios (~) and Mach Numbers (~)
A = [1,2,3, 12.8, 4.26];
M = 0:0.5:5;
%Define Aerodynamic Chord (C_mac,in) and Distance from Aerodynamic Chord (X_mac,in), %toc_mac is thickness over chord
c_mac = 13.3;
X_mac_Wing = 67.0;
toc_mac = 0.044;

XacOCmac_Surface_Matrix = zeros(length(M), length(A));
for i = 1:length(A)
    A1 = A(i);
    for j = 1:length(M)
        if M(j) >= 2
            XacOCmac_Surface = (A1 * (M(j)^2 - 1)^(1/2) - 0.67) / (2 * A1 * (M(j)^2 - 1)^(1/2) - 1);
        elseif M(j) <= 0.7
            XacOCmac_Surface = 0.25;
        end
        XacOCmac_Surface_Matrix(j, i) = XacOCmac_Surface;
    end
end

figure(1);
plot(M,XacOCmac_Surface_Matrix)
xlabel('Mach Number, M (~)');
ylabel('Surface Non-dimensional Aerodynamic Center,X_{AC}/C_{mac} (~)');
legend('A=1','A=2','A=3', 'GMLRS Tail Fin A=12.8', 'GMLRS Nose Fin A=4.26');
set(gca,"FontSize",20)
grid on;

ylim([0.2, 0.6]); 
hold off;

hold off;

%% Section 9
% % Hinge Moment Prediction
clear
clc
Cmac_w = 3.23; % Canard Mean Aerodynamic Center (in)
Cmac_t = 8.07; % Tail Mean Aerodynamic Center (in)
Aw = 4.26; % Canard aspect ratio
At = 12.79; % Tail Aspect ratio

M1 = 0.8;
M2 = 1.35;
M3 = 2;
M4 = 5;
q1 = 436; % dynamic pressure in psf at M = 0.8
q2 = 1242; % dynamic pressure in psf at M = 1.35
q3 = 2725; % dynamic pressure in psf at M = 2
q4 = 17031; % dynamic pressure in psf at M = 5
q = [q1 q2 q3 q4];
Mach = [M1 M2 M3 M4]; % Mach Numbers

AOA = 0:30; % alpha dash is the same as AOA
Ssurf_t = 833.97/144; %ft^2 tail planform area
Ssurf_w = 44.23/144; % ft^2 wing planform area
d = 8.94; % diameter in inches
Sref = (pi/4*d^2)/144; % Reference area in ft^2
S_w = Sref/Ssurf_w; % reference area / wing surface area
S_t = Sref/Ssurf_t; % reference area / tail surface area
CN_w = zeros(4,31);
CN_t = zeros(4,31);
N_w = zeros(4,31);
N_t = zeros(4,31);
HM_t = zeros(4,31);
HM_w = zeros(4,31);
for j = 1:31
    for i = 1:4
        if Mach(i) <= 1.35
            CN_w(i,j) = ((pi*Aw/2)*abs(sind(AOA(j))*cosd(AOA(j))))+2*((sind(AOA(j)))^2)*(1/S_w);
            CN_t(i,j) = ((pi*At/2)*abs(sind(AOA(j))*cosd(AOA(j))))+2*((sind(AOA(j)))^2)*(1/S_t);
            N_w(i,j) = CN_w(i,j)*(q(i)*Sref);
            N_t(i,j) = CN_t(i,j)*(q(i)*Sref);
            XAC_w = 0.25*Cmac_w;
            XAC_t = 0.25*Cmac_t;
            Xhl_t = 0.22*Cmac_t; % 22% of chord
            Xhl_w = 0.22*Cmac_w; % 22% of chord
            HM_w(i,j) = N_w(i,j)*(XAC_w-Xhl_w);
            HM_t(i,j) = N_t(i,j)*(XAC_t-Xhl_t);
            

        elseif Mach(i) == 2
            CN_w(i,j) = ((4*abs(sind(AOA(j))*cosd(AOA(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(AOA(j)))^2)*(1/S_w);
            CN_t(i,j) = ((4*abs(sind(AOA(j))*cosd(AOA(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(AOA(j)))^2)*(1/S_t);
            N_w(i,j) = CN_w(i,j)*(q(i)*Sref);
            N_t(i,j) = CN_t(i,j)*(q(i)*Sref);
            XAC_w = ((Aw*(sqrt(Mach(i)^2-1)-0.67))/(2*Aw*(sqrt(Mach(i)^2-1)-1)))*Cmac_w; 
            XAC_t = ((At*(sqrt(Mach(i)^2-1)-0.67))/(2*At*(sqrt(Mach(i)^2-1)-1)))*Cmac_t;
            Xhl_t = 0.25*Cmac_t;
            Xhl_w = 0.25*Cmac_w;
            HM_w(i,j) = N_w(i,j)*(XAC_w-Xhl_w);
            HM_t(i,j) = N_t(i,j)*(XAC_t-Xhl_t);
            
        elseif Mach(i) == 5
            CN_w(i,j) = ((4*abs(sind(AOA(j))*cosd(AOA(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(AOA(j)))^2)*(1/S_w);
            CN_t(i,j) = ((4*abs(sind(AOA(j))*cosd(AOA(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(AOA(j)))^2)*(1/S_t);
            N_w(i,j) = CN_w(i,j)*(q(i)*Sref);
            N_t(i,j) = CN_t(i,j)*(q(i)*Sref);
            XAC_w = ((Aw*(sqrt(Mach(i)^2-1)-0.67))/(2*Aw*(sqrt(Mach(i)^2-1)-1)))*Cmac_w;
            XAC_t = ((At*(sqrt(Mach(i)^2-1)-0.67))/(2*At*(sqrt(Mach(i)^2-1)-1)))*Cmac_t;
            Xhl_t = 0.25*Cmac_t;
            Xhl_w = 0.25*Cmac_w;
            HM_w(i,j) = N_w(i,j)*(XAC_w-Xhl_w);
            HM_t(i,j) = N_t(i,j)*(XAC_t-Xhl_t);
            
        end
    end
end

% % Plotting Figure 2.28
% Fig 2.28 data
Cmac1 = 13.3; % in
Xhl1 = 0.25*Cmac1; %in
Sref1 = 0.349; % ft^2
Sw1 = 2.55; %ft^2
h1 = 20*10^3; %ft
AOA1 = 0:30;
Sw = Sref1/Sw1;
Aw1 = 2.82; % baseline rocket AR

CN_w1 = zeros(4,31);
N_w1 = zeros(4,31);
HM_w1 = zeros(4,31);

for j = 1:31
    for i = 1:4
        if Mach(i) <= 1.35
            CN_w1(i,j) = ((pi*Aw1/2)*abs(sind(AOA(j))*cosd(AOA(j))))+2*((sind(AOA(j)))^2)*(1/Sw);          
            N_w1(i,j) = CN_w1(i,j)*(q(i)*Sref1);
            XAC_w1 = 0.25*Cmac1;            
            Xhl_w1 = 0.22*Cmac1; % 22% of chord
            HM_w1(i,j) = N_w1(i,j)*(XAC_w1-Xhl_w1);                  
        elseif Mach(i) == 2
            CN_w1(i,j) = ((4*abs(sind(AOA(j))*cosd(AOA(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(AOA(j)))^2)*(1/Sw);            
            N_w1(i,j) = CN_w1(i,j)*(q(i)*Sref1);
            XAC_w1 = ((Aw1*(sqrt(Mach(i)^2-1)-0.67))/(2*Aw1*(sqrt(Mach(i)^2-1)-1)))*Cmac1; 
            Xhl_w1 = 0.25*Cmac1;
            HM_w1(i,j) = N_w1(i,j)*(XAC_w1-Xhl_w1);           
        elseif Mach(i) == 5
            CN_w1(i,j) = ((4*abs(sind(AOA(j))*cosd(AOA(j))))/sqrt((Mach(i)^2)-1) + 2*(sind(AOA(j)))^2)*(1/Sw);            
            N_w1(i,j) = CN_w1(i,j)*(q(i)*Sref1);
            XAC_w1 = ((Aw1*(sqrt(Mach(i)^2-1)-0.67))/(2*Aw1*(sqrt(Mach(i)^2-1)-1)))*Cmac1; 
            Xhl_w1 = 0.25*Cmac1;
            HM_w1(i,j) = N_w1(i,j)*(XAC_w1-Xhl_w1);                      
        end
    end
end

% % Comparing Canard Hinge Moment with Figure 2.28
figure 


plot(AOA, HM_w(1,:), 'LineWidth',2)
hold on
plot(AOA, HM_w(2,:),'LineWidth',2)
hold on
plot(AOA, HM_w(3,:),'LineWidth',2)
hold on
plot(AOA, HM_w(4,:),'color', "#77AC30", 'LineWidth',2)
hold on

% Plotting Figure 2.28
plot(AOA, HM_w1(1,:),'-k','LineWidth',1)
hold on
plot(AOA, HM_w1(2,:), '--k', 'LineWidth',1)
hold on
plot(AOA, HM_w1(3,:),':k', 'LineWidth',2)
hold on
plot(AOA, HM_w1(4,:),'-.k', 'LineWidth',1)
hold on

lgd = legend ({'Mach 0.8 (GMLRS)', 'Mach 1.35 (GMLRS)', 'Mach 2 (GMLRS)', 'Mach 5 (GMLRS)', 'Mach 0.8 (Fig. 2.28)', 'Mach 1.35 (Fig. 2.28)', 'Mach 2 (Fig. 2.28)', 'Mach 5 (Fig. 2.28)' },'location','northeast', 'NumColumns',2);
fontsize(lgd,14,'points')
hold off
grid on
ax = gca;
ax.YAxis.Exponent = 0;
ylim([0 5000])
ylabel('Hinge Moment (HM) for canard, in-lb', 'FontSize',16)
xlabel("\alpha' = \alpha_{w} + \delta, Wing effective angle of attack, deg", 'FontSize',16)



% % Comparing Tail Hinge moment with figure 2.28
figure
plot(AOA, HM_t(1,:), 'LineWidth',2)
hold on
plot(AOA, HM_t(2,:),'LineWidth',2)
hold on
plot(AOA, HM_t(3,:),'LineWidth',2)
hold on
plot(AOA, HM_t(4,:),'color', "#77AC30", 'LineWidth',2)
hold on
% Plotting Figure 2.28
plot(AOA, HM_w1(1,:),'-k','LineWidth',1)
hold on
plot(AOA, HM_w1(2,:), '--k', 'LineWidth',1)
hold on
plot(AOA, HM_w1(3,:),':k', 'LineWidth',2)
hold on
plot(AOA, HM_w1(4,:),'-.k', 'LineWidth',1)
hold on
lgd = legend ({'Mach 0.8 (GMLRS)', 'Mach 1.35 (GMLRS)', 'Mach 2 (GMLRS)', 'Mach 5 (GMLRS)', 'Mach 0.8 (Fig. 2.28)', 'Mach 1.35 (Fig. 2.28)', 'Mach 2 (Fig. 2.28)', 'Mach 5 (Fig. 2.28)' },'location','northeast', 'NumColumns',2);
fontsize(lgd,14,'points')
hold off
grid on
ax = gca;
ax.YAxis.Exponent = 0;
ylim([0 30000])
ylabel('Hinge Moment (HM) for tail, in-lb', 'FontSize',16)
xlabel("\alpha' = \alpha_{w} + \delta, Wing effective angle of attack, deg", 'FontSize',16)
