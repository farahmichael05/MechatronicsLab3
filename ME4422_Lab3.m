clc
clear

%% Definitions 
%Everything with 1 is something I don't know the value of

%Lengths (in)
L1 = 1;
L2 = 9.45;
L3 = 1;
L4 = 1;
L5 = 1;
L6 = 0.75;
L7 = 0.45;
L8 = 3.7; %L2 in my notes (not in slack)
Lk_3 = 0.55;
Lk_4 = 1;


%Dimensions of platform
a_platform = 1.32;
b_platform = 1.25;

%Radius (in)
R = 0.05;

%Masses (UNITS)
M_platform = 1;
M_rack = 1;
Mk_3 = 1;
Mk_4 = 1;
M_E = 1; % Arms(on top of platform)
M_F = 1; % connecting arm with M_E

%Mass Moment of Inertia
J_pinion = 1; %Gear
J_2 = 1; %What is J_2?
J_leverarm = 1;

%Dampers
D_p1 = 1; %Friction with Pin 1
D_p2 = 1; %Friction with Pin 2
D_3 = 1;  %Friction with rack and ground
D_4 = 1;  %Friction with pinion and display
D_B = 1;  %Platform, M_E friction?

%Spring Stiffness

k_3 = 10.3*10^3;
k_4 = 5.04;
k_one_half = 1;

AppliedTorque = 1;

%% Equivalences

%Mass Equivalence -- Spring 2's location (Translational)

M_eq1 = M_platform + M_rack + (1/3)*Mk_3 + 2*(M_E + M_F); %add platform, spring 3, bars E and F (x2) and rack. THIS ASSUMES THEY ALL MOVE AT THE SAME VELOCITY (DOUBLE CHECK)
M_eq2 = J_pinion + (R^2)*M_eq1; %Lump to Pinion (Translation --> Rotation)

M_eq = (1/R^2)*(M_eq2) + (1/3)*(Mk_4); %Lump to k4 (Rotation --> Translation)

% Moment of Inertia Equivalence -- Pinion's location (Torsional)
%Comment consistency will be fixed later

J_eq1 = J_pinion + M_E*L1^2; % Lump bar E to pin between bar E and bar F
J_eq2 = J_2 + M_F*L4^2; % Lump bar F to ground
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

D_eq1R = D_B*L6^2 + D_p1; %Lump Platform damper to P, Translation to Rotation
D_eq2R = D_eq1R + D_p2; %add damper at S 
D_eq3T = (D_eq2R)* (1/L8^2); %to rack

D_eq = (D_eq3T)*(R^2) + D_4; %To pinion, Translation --> Rotation

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

% Plot angular displacement, velocity, and acceleration
figure;
grid on;
hold on;
title('Changes in Angular Position, Velocity, and Acceleration of Pinion over time')
fplot(solutionTheta, [0 10]);
fplot(solutionOmega, [0 10]);
fplot(solutionAlpha, [0 10]);
ylim([-2, 2]);
legend('\theta (rad)', '\omega (rad/s)', '\alpha (rad/s^2)');
xlabel('Time (s)')
hold off


%plotting x, v, a between times 0 and 10

figure;
grid on;
hold on;
title('Changes in Position, Velocity, and Acceleration of Pinion over time')
fplot(x, [0 10]);
fplot(v, [0, 10]);
fplot(a, [0 10]);
ylim([-0.06, 0.06]);
legend('x (in)', 'v (in/s)', 'a (in/s^2)');
xlabel('Time (s)')
hold off
