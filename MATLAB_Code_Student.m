%% ENM 056 : Machine Design Assignment
% V 1.0
% Responsible: ENM 056 Teaching group
% Contact: sharman@chalmers.se

% **************************************************************************
% This is a template file to help the studnets start up machine design 
% assignment in ENM 056. Please fill in the relevant sections one by one. 
% In order to run the code section by section, use "Run Section" option in
% MATLAB or press Ctrl + Enter. 
% Note: Please be careful when changing variable names in the template
% file. Also, the students may still need to write some sections by
% themselves.
% **************************************************************************

clear
close all
clc

%% Parameters of the reference machine
% Please fill the parameters of the machine below in SI units

mm = 1e-3; % mm to SI unit
OD_stator = 176 * mm; % Outer diameter of stator 
ID_stator = 124 * mm; % Inner diameter of stator
OD_rotor = 122 * mm ; % Outer diameter of rotor
ID_rotor = 60 * mm ; % Inner diameter of rotor
L_stack = 100 * mm  ; % Stack length
Hs0 = 0.5 * mm     ; % Slot opening height
Hs1 =  0.5 * mm     ; % Slot wedge height
Hs2 = 14 * mm      ; % Slot body height
w_tooth = 4.4 * mm  ; % Tooth width
Rs =   0.5 * mm     ; % Slot bottom radius fillet
Bs0 =  2 * mm     ; % Slot opening
N_pole =  8  ; % Number of poles
N_slot =  48  ; % Number of slots
w_layer = 2  ; % Number of winding layers
k_p =   5    ; % Coil pitch
D_wire = 0.75 * mm   ; % Wire diameter
f_Cu_max = 0.45 ; % Maximum Cu-fill factor
V_DC =  400    ; % DC link voltage [V]
mu_0 = 4 * pi * 1e-7; % Magnetic permeability of vacuum [H/m]
Omega_base = 4000;    % Base speed [rpm]
% Omega_max =         ; % Maximum speed [rpm]
N_parallel =  4      ; % Number of parallel branch
% N_strand =          ; % Number of strands
rho_Cu = 1.72e-8; % Resistivity of Cu [ohm/m]

%% Load line

% Start by assuming an air-gap flux density [T]
B_gap = 0.01: 0.001: 1;

% The effective air-gap cross section orthogonal to flux crossing the [sq.m]
% air-gap
A_gap = ; 

% Flux in the air-gap [Wb]
Phi_gap = ;

% Cross-section of tooth perpendicular to flux [sq.m]
A_tooth = ;

t_yoke = ; % yoke thickness [m]

% Cross-section of yoke perpendicul+ar to flux [sq.m]
A_yoke = ;

% Output values
disp('Cross-section of different parts')
fprintf('Air-gap cross-section = % .2f [mm^2] \n',A_gap * 1e6)
fprintf('Stator tooth cross-section = % .2f [mm^2] \n',A_tooth * 1e6)
fprintf('Stator yoke cross-section = % .2f [mm^2] \n \n',A_yoke * 1e6)

% Flux densities in different parts of the machine [T]
B_tooth = ; % Flux density in stator tooth
B_yoke = ; % Flux density in stator yoke 
B_rotor = ;

index = B_gap == 0.7;

% Output corresponding flux densities for air-gap flux density of 0.7 [T] in
% air-gap
fprintf('Flux densities in the different part of the machine when air-gap flux density is % .2f [T]\n', B_gap(index))
fprintf('Stator tooth flux density = % .2f [T] \n',B_tooth(index))
fprintf('Stator yoke flux density = % .2f [T] \n',B_yoke(index))
fprintf('Rotor yoke flux density = % .2f [T] \n \n',B_rotor(index))

% Import the B-H Curve of the M235-35A Steel from a TAB file
BH_data = importdata('SURA M235-35A - BH Curve @ 50 Hz.tab'); % import data
H_data = BH_data(:,1); % copy (Row All , Column One) as H
B_data = BH_data(:,2); % copy (Row All , Column Two) as B

% Calculated magnetic field intensity [A/m]
% Interpolation
method   = 'spline'; % 'linear' or 'spline' can be selected as interpolation method
H_rotor = interp1(B_data,H_data,B_rotor,method); % Interpolat rotor flux density, B_rotor to calculate corresponding H_rotor
H_stator_tooth = ; % Interpolat tooth flux density, B_tooth to calculate corresponding H_tooth
H_stator_yoke = ; % Interpolat tooth flux density, B_yoke to calculate corresponding H_yoke
H_gap = ; % Calculate magnetic field intensity in the air-gap. Note: Air-gap does not contain iron.

% Output magnetic field intensities
fprintf('Magnetic field intensities in the different part of the machine when air-gap flux density is % .2f [T] \n', B_gap(index))
fprintf('Air-gap Magnetic field intensity = % .2f [A/m] \n',H_gap(index))
fprintf('Stator tooth Magnetic field intensity = % .2f [A/m] \n',H_stator_tooth(index))
fprintf('Stator yoke Magnetic field intensity = % .2f [A/m] \n',H_stator_yoke(index))
fprintf('Rotor yoke Magnetic field intensity = % .2f [A/m] \n \n',H_rotor(index))

% Length of flux path [m]
l_rotor = ; % Length of flux path in rotor yoke
l_stator_yoke = ; % Length of flux path in stator yoke
l_stator_tooth = ; % Length of flux path in stator tooth
l_gap = ( ID_stator - OD_rotor ) / 2; % Length of flux path in air-gap

% Output values
disp('Length of flux path in different parts of the machine')
fprintf('Length of flux path in rotor = % .2f [mm] \n',l_rotor * 1e3)
fprintf('Length of flux path in stator yoke = % .2f [mm] \n',l_stator_yoke * 1e3)
fprintf('Length of flux path in stator tooth = % .2f [mm] \n',l_stator_tooth * 1e3)
fprintf('Length of flux path in air-gap = % .2f [mm] \n \n',l_gap * 1e3)

% MMF drop in different flux path [Aturn]
MMF_rotor = H_rotor * l_rotor; % MMF drop in rotor yoke
MMF_stator_yoke = ; % MMF drop in stator yoke
MMF_stator_tooth = ; % MMF drop in stator tooth
MMF_gap = ; % MMF drop in air-gap

% Output values
fprintf('MMF drop in different parts of the machine when air-gap flux density is % .2f [T] \n', B_gap(index))
fprintf('MMF drop in rotor yoke = % .2f [A-turn] \n',MMF_rotor(index))
fprintf('MMF drop in stator tooth = % .2f [A-turn] \n',MMF_stator_tooth(index))
fprintf('MMF drop in stator yoke = % .2f [A-turn] \n',MMF_stator_yoke(index))
fprintf('MMF drop in air-gap = % .2f [A-turn] \n\n',MMF_gap(index))

% Total MMF drop
MMF_total = ;

MMF_iron = ; % Total MMF drop in the iron parts
MMF_air = ; % Total MMF drop in air-gap

% Plot load line
figure(1)               % create Figure 1
clf                     % clear figure
subplot(1,2,1)
plot(MMF_total,Phi_gap*1e3, 'LineWidth', 2)
xlabel('Total MMF drop [A-turn]')
ylabel('Flux in air-gap[mWb]')
legend('Load line')
grid on

subplot(1,2,2)
hold on
plot(MMF_iron,Phi_gap*1e3, 'LineWidth', 2)
plot(MMF_air,Phi_gap*1e3, 'LineWidth', 2)
xlabel('MMF drop [A-turn]')
ylabel('Flux in air-gap[mWb]')
legend('MMF drop in iron','MMF drop in air')
grid on

%% Magnet dimension

t_mag = 1 * mm ; % Thickness of magnet segment
w_mag =  2 * mm; % Width of magnet segment

% Magnet data 
B_mag = []; 
H_mag = [];

MMF_mag = 2 * H_mag * t_mag; % MMF produced by magnet
Phi_mag = ; % Flux produced by the magnet

% Magnet characteristic
figure(2)               % create Figure 1
clf                     % clear figure
plot(MMF_mag,Phi_mag * 1e3, 'LineWidth', 2)
xlim([min(MMF_mag) 0])
xlabel('MMF produced by magnet [A-turn]')
ylabel('Flux due to magnet [mWb]')
legend('Demagnetization characteristic')
grid on

%% Finding the intersection
% To avoid extrapolation of data, we can limit the x-axis of the load line
% to maximum MMF that can be produced by magnet
index = MMF_total < max(abs(MMF_mag));
MMF_total_con = MMF_total(index);
Phi_gap_con = Phi_gap(index);

% Interpolate magnet characteristic corresponding to total MMF drop to have
% same x-axis. Don't forget the negative sine to move it to second quadrant
Phi_mag_interp = interp1( abs(MMF_mag), Phi_mag , MMF_total_con , method);
[value, index] = min(abs(Phi_mag_interp - Phi_gap_con));

Phi_gap_no_load = Phi_gap_con(index);
B_gap_no_load = Phi_gap_no_load / A_gap;
B_tooth_no_load = Phi_gap_no_load / A_tooth; % Flux density in stator tooth
B_yoke_no_load = Phi_gap_no_load / A_yoke; % Flux density in stator yoke

fprintf('Magnet thickness = % .2f [mm] \n',t_mag * 1e3)
fprintf('Magnet width = % .2f [mm] \n',w_mag * 1e3)
fprintf('No load flux in the air-gap = % .5f [Wb] \n',Phi_gap_no_load)
fprintf('No load flux density in the air-gap = % .2f [T] \n',B_gap_no_load)
fprintf('No load flux density in stator tooth = % .2f [T] \n',B_tooth_no_load)
fprintf('No load flux density in stator yoke = % .2f [T] \n\n',B_yoke_no_load)

%% Perform sensitivity analysis
% analysis 0: Do nothing
% analysis 1: Change in magnet thickness

analysis = 0;

switch analysis
    
    case 0
    % DO NOTHING    
    
    case 1       % Change in magnet thickness
        
        t_mag_1 = t_mag * 1.1; % 10% increase in magnet thickness
        t_mag_2 = t_mag * 0.9; % 10% decrease in magnet thickness
        
        MMF_mag_1 = 2 * H_mag * t_mag_1;
        MMF_mag_2 = 2 * H_mag * t_mag_2;
        
        % Magnet characteristic
        figure(4)               % create Figure 1
        clf                     % clear figure
        hold on
        plot(MMF_mag,Phi_mag * 1e3, 'LineWidth', 2)
        plot(MMF_mag_1,Phi_mag * 1e3, 'LineWidth', 2)
        plot(MMF_mag_2,Phi_mag * 1e3, 'LineWidth', 2)
        plot(-MMF_total,Phi_gap * 1e3, 'LineWidth', 2)
        hold off
        xlim([min(MMF_mag_1) 0])
        xlabel('MMF [A-turn]')
        ylabel('Flux [mWb]')
        legend('Magnet thickness 5 mm' , 'Magnet thickness 5.5 mm','Magnet thickness 4.5 mm', 'Load line')
        grid on
     
end

%% Slot area [sq.m]
% Tooth base [m]
Tooth_base = ;

% Tooth area [sq.m]
A_tooth = ;

% Tooth and slot area [sq.m]
A_tooth_slot = ;

% Slot area [sq.m]
A_slot = A_tooth_slot - A_tooth;

%% Stator winding design 

f_elec = ; % Electrical frequency [Hz]
omega_elec = ; % Electrical angular frequency [rad/s]
q = ; % Number of slots per pole per phase

% Total number of turns per phase
N_ph_total = ; 

% Number of series branch
N_series = ;

% Number of coils in series per phase
N_coil_ph = ;

% Number of turns per coil
N_turn_coil = round (N_ph_total / N_coil_ph);
N_ph_total = N_coil_ph * N_turn_coil;

% Number of conductors per slot
N_cond_slot = ;

% Cu area
A_Cu = ;

% Fill factor
f_Cu = A_Cu / A_slot;

% Maximum induced voltage at no-load
E_no_load = ;

% Maximum line to line induced voltage at no -load
E_LL_no_load = ;

% Maximum line to line induced voltage at max speed
E_LL_max_speed = ;

fprintf('Total number of turns per phase = % .0f \n',N_ph_total )
fprintf('Total number of coils in series per phase = % .0f \n',N_coil_ph )
fprintf('Total number of turns per coil = % .0f \n',N_turn_coil )
fprintf('Total number of conductors per slot = % .0f \n',N_cond_slot )
fprintf('Slot area = % .1f [mm^2] \n',A_slot * 1e6)
fprintf('Cu area = % .1f [mm^2] \n',A_Cu * 1e6)
fprintf('Fill factor = % .2f \n\n',f_Cu)

fprintf('Maximum per phase induced voltage at no-load = % .1f [V] \n',E_no_load)
fprintf('Maximum line to line induced voltage at no-load = % .1f [V] \n',E_LL_no_load)
fprintf('Maximum line to line induced voltage at %.0f [rpm] = % .1f [V] \n',Omega_max, E_LL_max_speed)
fprintf('Ratio of maximum line to line induced and DC link voltage at max speed of %.0f [rpm] = % .1f \n\n',Omega_max, E_LL_max_speed / V_DC)

%% Resistance and inductance calculation

R_phase = ;  % Resistance per phase [ohm]

Re_mag = ;
Re_gap = ;
Re_stator_tooth = ;
Re_stator_yoke = ;
Re_rotor_yoke = ;

Re_d = ;
Re_q = ;

N = ;

L_d = ;
L_q = ;

saliency = L_q / L_d;

fprintf('Resistance per phase = % .1f [mOhm] \n',R_phase*1e3)
fprintf('Reluctance of d-flux path = % .1f [H^-1] \n',Re_d)
fprintf('Reluctance of q-flux path = % .1f [H^-1] \n',Re_q)
fprintf('D-axis inductance = % .1f [mH] \n',L_d*1e3)
fprintf('Q-axis inductance = % .1f [mH] \n',L_q*1e3)
fprintf('Saliency ratio = % .1f \n',saliency)
