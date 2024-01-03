%% Section 10
clc
close 
clear 

% Code Notes
% Code to solve for the AC loaction as a ratio of the distance to the body
% aerodynamic cenetr to the length of the nose of the missile

% % Section 1

% Variables

S_Surface = 1; % Surface Area (ft^2)
S_Ref = 0.435; % Reference Area (ft^2)
n_Surface = 8; % Number of Surface Planforms (~)
q = 2725; % Dynamic Pressure (lb/ft^2)
c_mac = 1; % Length of Mean Areo Chord (ft)
M = 0:0.1:6; % Mach Number (~)
MQC = [0.00001 0.0001 0.001 0.01]; % (ft/lb)
MQCGMLRS = M./(q*c_mac); %

% Equations

NSS = 0:0.1:50; %
NSSGMLRS = (n_Surface*S_Surface)/(S_Ref); %
C_D_0_SurfaceGMLRS = (0.0133.*(MQCGMLRS.^0.2)).*2.*(NSSGMLRS); %
C_D_0_Surface1 = (0.0133.*(MQC(1,1).^0.2)).*2.*(NSS); %
C_D_0_Surface2 = (0.0133.*(MQC(1,2).^0.2)).*2.*(NSS);
C_D_0_Surface3 = (0.0133.*(MQC(1,3).^0.2)).*2.*(NSS);
C_D_0_Surface4 = (0.0133.*(MQC(1,4).^0.2)).*2.*(NSS);

% Graph
figure
%plot(NSSGMLRS,C_D_0_SurfaceGMLRS,'k')
hold on
plot(NSS,C_D_0_Surface1,'b')
plot(NSS,C_D_0_Surface2,'r')
plot(NSS,C_D_0_Surface3,'m')
plot(NSS,C_D_0_Surface4,'g')
hold off
ylim([0 0.4])
xlim([0 50])
yticks([0 0.1 0.2 0.3 0.4]);
xticks([0 10 20 30 40 50]);
ylabel('(C_{D_{0}})_{Surface,Friction} (~)');
xlabel('nS_{Surface}/S_{Ref} (~)');
legend('M/(q cmac) = GMLRS ft/lb','M/(q cmac) = 0.00001 ft/lb','M/(q cmac) = 0.0001 ft/lb','M/(q cmac) = 0.001 ft/lb','M/(q cmac) = 0.01 ft/lb','','Location','northwest');


% % Section 2

% Variables

n_W = 4; % Number of Wings (~)
n_WC = 4; % Number of Canards (~)
n_WT = 4; % Number of Tail Wings (~)
delta_LE = 2; % Leading Edge Section Angle (deg)
delta_LEC = 2; % Canard Leading Edge Section Angle (deg)
delta_LET = 2; % Tail Leading Edge Section Angle (deg)
gamma_LE = [5 10 20 90]; % Canard Leading Edge Sweep Angle (deg)
gamma_LEC = 87; % Canard Leading Edge Sweep Angle (deg)
gamma_LET = 90; % Tail Leading Edge Sweep Angle (deg)
t_mac = 0.0208; % Max Thickness of MAC (ft)
t_macC = 0.0208; % Canard Max Thickness of MAC (ft)
t_macT = 0.0417; % Tail Max Thickness of MAC (ft)
b = 0.269; % Span (ft)
bC = 0.269; % Canard Span (ft)
bT = 0.585; % Tail Span (ft)
S_Ref = 0.435; % Reference Area (ft^2)
%M = 1:0.1:6; % Mach Number (~)

% Equations

M_gamma_LE1 = M.*cosd(gamma_LE(1,1)); % Mach Number Perpendicular to
% Leading Edge (~)
M_gamma_LE2 = M.*cosd(gamma_LE(1,2)); % Mach Number Perpendicular to
% Leading Edge (~)
M_gamma_LE3 = M.*cosd(gamma_LE(1,3)); % Mach Number Perpendicular to
% Leading Edge (~)
M_gamma_LE4 = M.*cosd(gamma_LE(1,4)); % Mach Number Perpendicular to
% Leading Edge (~)
M_gamma_LEC = M.*cosd(gamma_LEC); % Mach Number Perpendicular to Leading
% Edge of Canard (~)
M_gamma_LET = M.*cosd(gamma_LET); % Mach Number Perpendicular to Leading
% Edge of Tail (~)

C_D_0_Surface_Wave1 = n_W.*(1.429./(M_gamma_LE1.^2)).*(((1.2.* ...
    (M_gamma_LE1).^2).^3.5).*((2.4./(2.8.*(M_gamma_LE1.^2)-0.4)) ...
    .^2.5)-1).*(((sind(delta_LE).^2).*cosd(gamma_LE(1,1)).*t_mac.* ...
    (b./(S_Ref))));
% Wing Surface Wave Zero Lift Drag Coefficient (%)

C_D_0_Surface_Wave2 = n_W.*(1.429./(M_gamma_LE2.^2)).*(((1.2.* ...
    (M_gamma_LE2).^2).^3.5).*((2.4./(2.8.*(M_gamma_LE2.^2)-0.4)) ...
    .^2.5)-1).*(((sind(delta_LE).^2).*cosd(gamma_LE(1,2)).*t_mac.* ...
    (b./(S_Ref))));
% Wing Surface Wave Zero Lift Drag Coefficient (%)

C_D_0_Surface_Wave3 = n_W.*(1.429./(M_gamma_LE3.^2)).*(((1.2.* ...
    (M_gamma_LE3).^2).^3.5).*((2.4./(2.8.*(M_gamma_LE3.^2)-0.4)) ...
    .^2.5)-1).*(((sind(delta_LE).^2).*cosd(gamma_LE(1,3)).*t_mac.* ...
    (b./(S_Ref))));
% Wing Surface Wave Zero Lift Drag Coefficient (%)

C_D_0_Surface_Wave4 = n_W.*(1.429./(M_gamma_LE4.^2)).*(((1.2.* ...
    (M_gamma_LE4).^2).^3.5).*((2.4./(2.8.*(M_gamma_LE4.^2)-0.4)) ...
    .^2.5)-1).*(((sind(delta_LE).^2).*cosd(gamma_LE(1,4)).*t_mac.* ...
    (b./(S_Ref))));
% Wing Surface Wave Zero Lift Drag Coefficient (%)

C_D_0_Surface_WaveC = n_WC.*(1.429./(M_gamma_LEC.^2)).*(((1.2.* ...
    (M_gamma_LEC).^2).^3.5).*((2.4./(2.8.*(M_gamma_LEC.^2)-0.4)) ...
    .^2.5)-1).*(((sind(delta_LEC).^2).*cosd(gamma_LEC).*t_macC.* ...
    (bC./(S_Ref))));
% Canard Surface Wave Zero Lift Drag Coefficient (%)

C_D_0_Surface_WaveT = n_WT.*(1.429./(M_gamma_LET.^2)).*(((1.2.* ...
    (M_gamma_LET).^2).^3.5).*((2.4./(2.8.*(M_gamma_LET.^2)-0.4)) ...
    .^2.5)-1).*(((sind(delta_LET).^2).*cosd(gamma_LET).*t_macT.* ...
    (bT./(S_Ref))));
% Tail Surface Wave Zero Lift Drag Coefficient (%)

AAC = (C_D_0_Surface_WaveC./(n_WC.*cosd(gamma_LEC).*t_macC*(bC/S_Ref)))+C_D_0_Surface1;
AAT = (C_D_0_Surface_WaveT./(n_WT.*cosd(gamma_LET).*t_macT*(bT/S_Ref)))+C_D_0_Surface1;
AA1 = (C_D_0_Surface_Wave1./(n_W.*cosd(gamma_LE(1,1)).*t_mac*(b/S_Ref)))+C_D_0_Surface1;
AA2 = (C_D_0_Surface_Wave2./(n_W.*cosd(gamma_LE(1,2)).*t_mac*(b/S_Ref)))+C_D_0_Surface1;
AA3 = (C_D_0_Surface_Wave3./(n_W.*cosd(gamma_LE(1,3)).*t_mac*(b/S_Ref)))+C_D_0_Surface1;
AA4 = (C_D_0_Surface_Wave4./(n_W.*cosd(gamma_LE(1,4)).*t_mac*(b/S_Ref)))+C_D_0_Surface1;


% Graph

figure
semilogy(M_gamma_LEC,abs(AAC),'b')
hold on
semilogy(M_gamma_LET,abs(AAT),'k')
semilogy(M_gamma_LE1,abs(AA1),'m')
semilogy(M_gamma_LE2,abs(AA2),'g')
semilogy(M_gamma_LE3,abs(AA3),'r')
semilogy(M_gamma_LE3,abs(AA4),'c')
hold off
ylim([0.01 10])
xlim([0 6])
yticks([0.01 0.1 1 10]);
xticks([0 1 2 3 4 5 6]);
ylabel('(C_{D_{0}})_{Surface,Friction}');
xlabel('nS_{Surface}/S_{Ref}');
legend('Canard','Tail','DeltaLE = 5{\circ}','DeltaLE = 10{\circ}','DeltaLE = 20{\circ}','DeltaLE = 90{\circ} (Blunt)','Location','northwest');


% % Section 3

% Variables

alphai = 0:0.1:60; % Wing Local Angle of Attack (deg)
C_D_0 = [0.01 0.005 0.005]; % Zero-Lift Drag Coefficient (~)
C_D_0C = 0.01; % Canard Zero-Lift Drag Coefficient (~)
C_D_0T = 0.01; % Tail Zero-Lift Drag Coefficient (~)
A = [2 2 4]; % Apsect Ratio (~)
AC = 2.7; % Canard Apsect Ratio (~)
AT = 3.2; % Tail Apsect Ratio (~)
%bC = 1; % Span (in)
%bT = 1; % Span (in)
%S_Ref = 0.435; % Reference Area (ft^2)

% Equations

C_N_1 = (((pi*A(1,1))/2).*abs(sind(alphai).*cosd(alphai)))+(2.*(sind(alphai).^2)); %
LD1 = ((C_N_1.*cosd(alphai))-(C_D_0(1,1).*sind(alphai)))./((C_N_1.*sind(alphai))+(C_D_0(1,1).*cosd(alphai))); %

C_N_2 = (((pi*A(1,2))/2).*abs(sind(alphai).*cosd(alphai)))+(2.*(sind(alphai).^2)); %
LD2 = ((C_N_2.*cosd(alphai))-(C_D_0(1,2).*sind(alphai)))./((C_N_2.*sind(alphai))+(C_D_0(1,2).*cosd(alphai))); %

C_N_3 = (((pi*A(1,3))/2).*abs(sind(alphai).*cosd(alphai)))+(2.*(sind(alphai).^2)); %
LD3 = ((C_N_3.*cosd(alphai))-(C_D_0(1,3).*sind(alphai)))./((C_N_3.*sind(alphai))+(C_D_0(1,3).*cosd(alphai))); %

C_N_C = (((pi*AC)/2).*abs(sind(alphai).*cosd(alphai)))+(2.*(sind(alphai).^2)); %
LDC = ((C_N_C.*cosd(alphai))-(C_D_0C.*sind(alphai)))./((C_N_C.*sind(alphai))+(C_D_0C.*cosd(alphai))); %

C_N_T = (((pi*AT)/2).*abs(sind(alphai).*cosd(alphai)))+(2.*(sind(alphai).^2)); %
LDT = ((C_N_T.*cosd(alphai))-(C_D_0T.*sind(alphai)))./((C_N_T.*sind(alphai))+(C_D_0T.*cosd(alphai))); %


% Graph
figure
plot(alphai,LDC,'b')
hold on
plot(alphai,LDT,'k')
plot(alphai,LD1,'m')
plot(alphai,LD2,'g')
plot(alphai,LD3,'r')
hold off
ylim([0 20])
xlim([0 60])
yticks([0 4 8 12 16 20]);
xticks([0 20 40 60]);
ylabel('Canard Lift-to-Drag Ratio, L/D (~)');
xlabel('\alpha`, Canard Local Angle of Attack (\circ)');
legend('GMLRS Canard','GMLRS Tail','C_{D_{0}}=0.01, A=2','C_{D_{0}}=0.005, A=2','C_{D_{0}}=0.005, A=4','Location','northeast');

%% Section 11
clc
clear

Mcritical = 0:0.1:2; %Critical mach number (~)
AstarOA = (216/125).*Mcritical.*(1+Mcritical.^2./5).^(-3); %choked flow area over area
toh = (1-AstarOA)./2; %thickness over height

tNose_GMLRS = 0.06;%in, thickness of nose fin
hNose_GMLRS = 1;%in, height of nose fin
tohNose_GMLRS = tNose_GMLRS/hNose_GMLRS; %thickness/height of nose fin
tTail_GMLRS = 0.06;%in, thickness of tail fin
hTail_GMLRS = 2.5;%in, height of tail fin
tohTail_GMLRS = tTail_GMLRS/hTail_GMLRS; %thickness/height of tail fin

AstarONose_GMLRS = 1-2*(tohNose_GMLRS);
AstarOTail_GMLRS = 1-2*(tohTail_GMLRS);

figure;
plot(toh,Mcritical)
hold on
xline(tohNose_GMLRS, 'r')
hold
xline(tohTail_GMLRS,'g')
hold on
xlim([0 0.1])
ylim([0 2])
legend('Lattice Fin Choked Flow','Nose Fin','Tail Fin', 'FontSize', 30)
xlabel('Cell Wall Thickness/Cell Height, t/h (~)', 'FontSize', 30)
ylabel('Mach Number for Choked Flow, M_{Critical} (~)', 'FontSize', 30)
set(gca, 'FontSize', 30)

%% Section 12
clc
clear
%.........................Section 12................................
AR_C = 2.105;% Canard Aspect Ratio
AR_T = 3.33;%Tail Fin Aspect Ratio
d = 8.94;%Diameter (Inches)
M = 0:0.05:5;%Mach number
S_C = 15.7;%Canard Area of a pair (in2)
S_T = 73.3;%Tail wet Area of a pair(in2)
S_ref = 62.77;% Reference Area (in2)
X_AC_B = 12.28;%Location of Aerodynamic Center (inches)
X_CG_Full = 78.4;%Location of Center of Gravity @Full Fuel Load (inches)
X_CG_Spent = 72.5;%Location of Center of Gravity @Full Fuel Load (inches)
MAC_C = 2.87;%Length of Mean Aerodynamic Chord of Canard(inches)
MAC_T = 7.54;%Length of Mean Aerodynamic Chord of Tail Fin(inches)
X_MAC_C = 7.18;%Location of Canard's Mean Aerodynamic Chord(inches) 
X_MAC_T = 138.67;%Location of Canard's Mean Aerodynamic Tail Fin(inches)
%........................AC of Canard and Tail Fin..................
M = 0:1:5;% Mach Number
A = size(M);
X_AC_C = zeros(1,A(:,2));
X_AC_T = zeros(1,A(:,2));
%X_AC for Canard section
for i = 1:length(M)
    if M(i) <= 0.7
        X_AC_C = 0.25;
    elseif M(i) >= 2
            X_AC_C = ((sqrt((M.^2)-1).*AR_C)-0.67)/((sqrt((M.^2)-1)*2*AR_C)-1);
    end
end
    X_C_AC(1,i) = X_AC_C;

%X_AC for Tail fin section
for i = 1:length(M)
    if M(i) <= 0.7
        X_AC_T = 0.25;
    else if M(i) >= 2
            X_AC_T = ((sqrt((M.^2)-1).*AR_T)-0.67)/((sqrt((M.^2)-1)*2*AR_T)-1);
    end
    end
    X_T_AC(1,i) = X_AC_T;
end
figure(1)
plot(M,X_C_AC)
xlabel('Mach Number, M (~)')
ylabel({'Canard Aerodynamic Center Location over Canard Mean Dynamic Chord,'; '(X_AC)canard/(MAC)canard (~)'})
ylim([0.2, 0.6]);
set(gca,'FontSize',18)
set(gcf,'color','w')
figure(2)
plot(M,X_T_AC)
ylim([0.2, 0.6]);
xlabel('Mach Number, M (~)')
ylabel({'Tail Fin Aerodynamic Center Location over Tail Fin Mean Dynamic Chord,'; '(X_AC)Fin/(MAC)Fin (~)'})
set(gca,'FontSize',18)
set(gcf,'color','w')
%........................final Xmac/Cmac for Canard............
M = 0:0.05:5;%Mach Number
X_AC_1 = zeros(1,21);
X_AC_1(1,:) = 0.25;
M_2_5 = 2.05:0.05:5;% Mach Number from 2 to 5
A = size(M_2_5);
X_AC_2 = zeros(1,A(:,2));
for i = 1:length(M_2_5)
    X_AC = ((sqrt((M_2_5(i)^2)-1)*AR_C)-0.67)/((sqrt(((M_2_5(i))^2)-1)*2*AR_C)-1);
    X_AC_2(1,i) = X_AC;
end
M_1_2 = 1.05:0.05:2;% Mach Number from 1.05 to 2
A = size(M_1_2);
X_AC_1_2 = zeros(1,A(:,2));
for i = 1:length(M_1_2)
        X_AC = 0.25+0.0119658*(i-1); 
    X_AC_1_2(1,i) = X_AC;
end
X_AC_C = [X_AC_1(1,:),X_AC_1_2(1,:),X_AC_2(1,:)];
figure(3)
plot(M,X_AC_C)
ylim([0.2, 0.6]);
xlabel('Mach Number, M (~)')
ylabel({'Canard Aerodynamic Center Location over Mean Dynamic Chord,'; '(X_AC)canard/(MAC)canard (~)'})
ylim([0.2, 0.6]);
set(gca,'FontSize',14)
set(gcf,'color','w')
%........................final Xmac/Cmac for Tail Fin............
X_AC_1 = zeros(1,21);
X_AC_1(1,:) = 0.25;
M_2_5 = 2.05:0.05:5;% Mach Number from 2 to 5
A = size(M_2_5);
X_AC_2 = zeros(1,A(:,2));
for i = 1:length(M_2_5)
    X_AC = ((sqrt((M_2_5(i)^2)-1)*AR_T)-0.67)/((sqrt(((M_2_5(i))^2)-1)*2*AR_T)-1);
    X_AC_2(1,i) = X_AC;
end
M_1_2 = 1.05:0.05:2;% Mach Number from 1.05 to 2
A = size(M_1_2);
X_AC_1_2 = zeros(1,A(:,2));
for i = 1:length(M_1_2)
        X_AC = 0.25+0.01216335*(i-1); 
    X_AC_1_2(1,i) = X_AC;
end
X_AC_T = [X_AC_1(1,:),X_AC_1_2(1,:),X_AC_2(1,:)];
figure(4)
plot(M,X_AC_T)
ylim([0.2, 0.6]);
xlabel('Mach Number, M (~)')
ylabel({'Tail Fin Aerodynamic Center Location over Mean Dynamic Chord,'; '(X_AC)Fin/(MAC)Fin (~)'})
set(gca,'FontSize',14)
set(gcf,'color','w')
%..................................CNa.................................
C_N_AOA_B = 2;%  %body's Normal Force Coefficient Derivative From Angle of Attack (/rad) 
M_Critiical_C = sqrt(1+(8/(pi*AR_C))^2);
M_Critiical_T = sqrt(1+(8/(pi*AR_T))^2);
A = size(M);
C_N_AoA_C = zeros(1,A(:,2));
C_N_AoA_T = zeros(1,A(:,2));
for i = 1:length(M)
    if M(i) > M_Critiical_C
        C_N_AoA = (4/sqrt((M(i)^2)-1))*(S_C/S_ref);
    else if M(i) < M_Critiical_C
            C_N_AoA = (pi*AR_C/2)*(S_C/S_ref);
    end
    end
    C_N_AoA_C(1,i) = C_N_AoA;
end
for i = 1:length(M)
    if M(i) > M_Critiical_T
        C_N_AoA = (4/sqrt((M(i)^2)-1))*(S_T/S_ref);
    else if M(i) < M_Critiical_T
            C_N_AoA = (pi*AR_C/2)*(S_T/S_ref);
    end
    end
    C_N_AoA_T(1,i) = C_N_AoA;
end
%.......................Calculation for the required tail area.....
X_AC_T = X_AC_T*MAC_T+X_MAC_T;%location of aerodynamic center of tail fin (inches)
X_AC_C = X_AC_C*MAC_C+X_MAC_C;%location of aerodynamic center of Canard (inches)
b_Full = C_N_AOA_B*(X_CG_Full-X_AC_B)/d;% (XCG-(XAC)B)/d at full fuel load(~)
b_Spent = C_N_AOA_B*(X_CG_Spent-X_AC_B)/d;% (XCG-(XAC)B)/d at spent fuel load(~)
C_Full = C_N_AoA_C*(S_C/S_ref).*(X_CG_Full-X_AC_C)/d;% (XCG-(XAC)C)/d at full fuel load(~)
C_Spent = C_N_AoA_C*(S_C/S_ref).*(X_CG_Spent-X_AC_C)/d;% (XCG-(XAC)C)/d at spent fuel load(~)
T_Full = (d*S_ref)./((X_AC_T-X_CG_Full).*C_N_AoA_T);% (d/((XAC)T-XCG) at full fuel load(~)
T_Spent = (d*S_ref)./((X_AC_T-X_CG_Spent).*C_N_AoA_T);% (d/((XAC)T-XCG) at Spent fuel load(~)
ST_Full = (b_Full+C_Full).*T_Full;
ST_Spent = (b_Spent+C_Spent).*T_Spent;
figure(5)
ST1 = plot(M,ST_Spent/S_ref,'k');
hold on
ST2 = plot(M,ST_Full/S_ref,'r');
hold off
set(gca,'FontSize',20)
set(gcf,'color','w')
xlabel('Mach Number, M (~)')
ylabel({'tability Tail Area/Reference Area';'(ST)netural/(Sref) (~)'})
ST3 = legend([ST1 ST2],{'GMLRS Rocket At Spent Fuel Load','GMLRS Rocket AT Full Fuel Load'},'Location','northwest');

%% Section 13

clear
clc
Ac = 2.1; % Canard aspect ratio
At = 3.33; % Tail Aspect ratio
M = 2.5;
AOA1 = 5; % lower AOA 
AOA2 = 20; % higher AOA
Ssurf_t = 5.79; %ft^2 tail wet area
Ssurf_c = 73.3/144; % ft^2 canard planform area
d = 8.94; % diameter in inches
Sref = (pi/4*d^2)/144; % Reference area in ft^2
S_w = Sref/Ssurf_c; % reference area / canard surface area
S_t = Sref/Ssurf_t; % reference area / tail surface area

% for AOA 1 
CN_c1 = ((4*abs(sind(AOA1)*cosd(AOA1)))/sqrt((M^2)-1) + 2*(sind(AOA1))^2)*(1/S_w);            
CN_t1 = ((4*abs(sind(AOA1)*cosd(AOA1)))/sqrt((M^2)-1) + 2*(sind(AOA1))^2)*(1/S_t);            
CN_body1 = 2*AOA1/57.3; 


% for AOA 2
CN_c2 = ((4*abs(sind(AOA2)*cosd(AOA2)))/sqrt((M^2)-1) + 2*(sind(AOA2))^2)*(1/S_w);            
CN_t2 = ((4*abs(sind(AOA2)*cosd(AOA2)))/sqrt((M^2)-1) + 2*(sind(AOA2))^2)*(1/S_t); 
CN_body2 = 2*AOA2/57.3;         

% Total Normal force coefficient 
CN_total1 = CN_c1+CN_t1+CN_body1; % For AOA = 5
CN_total2 = CN_c2+CN_t2+CN_body2; % For AOA = 20

%Plotting bar graph for the Normal force Coefficient distribution 
x = [5, 20];
y = [CN_body1 CN_t1 CN_c1; CN_body2 CN_t2 CN_c2];
bar(x,y,'stacked')
legend('(C_N)_{Body}','(C_N)_{Tail}','(C_N)_{Canard}','FontSize',18);
xlabel('Angle of Attack, \alpha, deg','FontSize',18)
ylabel('Normal Force Coefficient of GMLRS, C_N (~)','FontSize',18)


