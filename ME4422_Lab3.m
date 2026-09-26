clc
clear

%% Definitions 
%Everything with 1 is something I don't know the value of

%Lengths (cm)
L1 = 1;
L2 = 24.003;
L3 = 1;
L4 = 1;
L5 = 1;
L6 = 1.905;
L7 = 2.286;
L8 = 9.398; %L2 in my notes (not in slack)
Lk_3 = 1.397;
Lk_4 = 1;

%Dimensions of platform
a_platform = 3.3528;
b_platform = 3.175;

%Radius (cm)
R = 0.127;

%Masses (grams)
M_platform = 9.7902;
M_rack = 7.60;
Mk_3 = 7.35;
Mk_4 = 0.29;
M_E = 0; % negligible
M_F = 0; % negligible


%Mass Moment of Inertia (g*cm^2)
J_pinion = 0.0021017; %Pinion
J_DC = 357.17; %Bar DC
J_leverarm = 8.2429;
J_AB = 2645.84; %bar AB

%Dampers
D_P = 0.3016; %Friction with Pin P
D_R = 0;  %Friction with rack, negligible 
D_pin = 0; %Friction with pinion, negligible
D_D = 9.425;  %friction at pivot D
D_A = 9.425;  %friction at pivot A
D_B = 9.425;  %Platform

%Spring Stiffness

k_3 = 103;
k_4 = 0.0504 ;
k_one_half = 1;

AppliedTorque = 0.1;

%% Equivalences

%Mass Equivalence -- Spring 2's location (Translational)

M_eq1 = (M_platform + (1/3)*Mk_3); %Add mass of platform and spring 3
JM_eq1 = (M_eq1)*L2^2 + 2*(J_AB); %Lump Meq1 to bar AB (Translation --> Rotation)
JM_eq2 = JM_eq1*(L3^2)/(L2^2) + 2*J_DC; %lump JMeq1 to bar DC
JM_eq3 = ((L6/L2)^2)*JM_eq2 + J_leverarm; %lump JMeq2 to the leverarm
M_eq2 = (JM_eq3)/L6^2 + M_rack; %lump JMeq3 to the rack
JM_eq4 = (M_eq2*R^2) + J_pinion; %add the pinion (Translation --> Rotation)

M_eq = (JM_eq4)/R^2 + (1/3)*(Mk_4); %Lump to k4 (Rotation --> Translation)

% Moment of Inertia Equivalence -- Pinion's location (Torsional)

J_eq1 = J_AB + M_E*L1^2; % Lump bar E to pin between bar E and bar F
J_eq2 = J_DC + M_F*L4^2; % Lump bar F to ground
MJ_eq1 = J_eq2 * 1/(L3^2); % Lump Jeq2 to C
J_eq3 = (L5^2)*MJ_eq1; % C to A
MJ_eq2 = (J_eq1 + J_eq3)/(L2^2); %Add Jeq1 and Jeq3 and lump to B
MJ_eq3 = 2*MJ_eq2; % both sides
MJ_eq4 = MJ_eq3 + M_platform; 
J_eq4 = (L6)^2*(MJ_eq4) + J_leverarm; %Lump B to P
MJ_eq5 = J_eq4/(L7^2) + M_rack; %lump P to S and add mass of rack

J_eq = (R^2)*(MJ_eq5) + J_pinion; %lump to pinion

%Spring Stiffness Equivalence -- Pinion's Location (Torsional)

k_eq1R = (k_one_half) * L1^2; %Lump k/2 to A
k_eq2R = (k_one_half) * L4^2; %Lump k/2 (second) to F to D
k_eq3T = k_eq2R/(L3^2); % Lump k_eq2R from D to C
k_eq4R = (L5^2)*k_eq3T; %k_eq3T from C to A
k_eq5T = (k_eq4R + k_eq1R)/(L2^2); % add k_eq1R and k_eq4R (make translation)
k_eq6T = 2*k_eq5T; % Both sides
k_eq7T = k_eq6T + k_3; % add k_3
k_eq8R = (L6^2)*k_eq7T; %move k_eq7T from B to P
k_eq9T = k_eq8R/(L7^2) + k_4; %move to S and add k_4

k_eq = (R^2)*k_eq9T; %move to pinion

%Damper Equivalence -- Pinion's Location (Torsional)

D_eq1 = D_D/(L3^2); %Lump D to C
D_eq2 = D_eq1 * (L5^2); %Lump C to A
D_eq3 = (D_A + D_eq2)/(L2^2); %Add A and move to B
D_eq4 = 2*D_eq3; %Both Sides
D_eq5 = D_eq4 + D_B;%Add platform damper
D_eq6 = D_eq5*(L6^2) + D_P; %Lump B to P
D_eq7 = D_eq6/(L7^2) + D_R; %Lump P to S and add Rack

D_eq = D_eq7*(R^2) + D_pin; %Lump to Pinion

%% Model Solver and Plotting

syms theta(t)

eqn = J_eq*diff(theta, t, 2) + D_eq*diff(theta, t, 1) +k_eq*theta == AppliedTorque;
Dtheta = diff(theta, t);

%Initial Conditions
initialCond = [theta(0) == 0, Dtheta(0) == 0];

%Solving for theta(t)
solutionTheta = dsolve(eqn, initialCond);
solutionOmega = diff(solutionTheta);
solutionAlpha = diff(solutionOmega);

%finding x v and a

x = R*solutionTheta;
v = R*solutionOmega;
a = R*solutionAlpha;

% Angular quantities figure (subplots)
figure;
sgtitle('Changes in Angular Position, Velocity, and Acceleration of Pinion over time')

subplot(3,1,1)
fplot(solutionTheta, [0 20],'r', 'LineWidth', 1.2);
ylim([-0.001, 0.2]);
ylabel('\theta (rad)')
grid on

subplot(3,1,2)
fplot(solutionOmega, [0 20],'b', 'LineWidth', 1.2);
ylim([-0.2, 0.2]);
ylabel('\omega (rad/s)')
grid on

subplot(3,1,3)
fplot(solutionAlpha, [0 20], 'c', 'LineWidth', 1.2);
ylim([-0.35, 0.35]);
ylabel('\alpha (rad/s^2)')
xlabel('Time (s)')
grid on

% Translational quantities figure (subplots)
figure;
sgtitle('Changes in Position, Velocity, and Acceleration of Rack over time')

subplot(3,1,1)
fplot(x, [0 20],'r', 'LineWidth', 1.2);
ylim([-0.0001, 0.03])
ylabel('x (cm)')
grid on

subplot(3,1,2)
fplot(v, [0 20],'b', 'LineWidth', 1.2);
ylim([-0.03, 0.03]);
ylabel('v (cm/s)')
grid on

subplot(3,1,3)
fplot(a, [0 20],'c', 'LineWidth', 1.2);
ylim([-0.1, 0.1]);
ylabel('a (cm/s^2)')
xlabel('Time (s)')
grid on
