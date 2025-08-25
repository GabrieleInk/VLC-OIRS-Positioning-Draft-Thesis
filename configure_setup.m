function configure_setup

clear variables;

% Parameters setting
Psi     = 70;               % LED half-power semiangle [degree]
rho     = 0.95;             % Reflection coefficient
%A_pd    = 1e-04;           % 1cm^2 - Physical area of the PD [m^2] [C2022]
A_pd    = 0.2e-04;          % Physical area of the PD [m^2]
T_of    = 1;                % Optical Filter Gain
a       = 1.5;              % Refractive index
Phi_FoV = 70;               % Field of view [degree]
B       = 5e6;              % System bandwidth [Hz]
%R_pd    = 2.2e-8;          % Responsivity [A][m^2] / [W] -> https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6868970
R_pd    = 0.54;             % (Sensitivity) [A/W]
q       = 3;                % Conversion ration of optical-to-electrical power
N       = 10e-21;           % Power spectral density [A^2/Hz]
lumen_level = 1000;
p = lumen_level / 683;      % Transmission power 6000 [Lumens] -> [Watt]

% Noise parameter
q_0     = 1.602e-19;        % Electronic charge [Coulombs]
I_bg    = 5e-12;            % Background light current 5 [pA]
k_B     = 1.38064852e-23;   % Boltzmann constant [Joule/Kelvin]
T_k     = 295;              % Absolute temperature [K]
G_0     = 10;               % Open-loop voltage gain
eta     = 1.12e-6;          % Fixed capacitance of photo detector per unit area [F/m^2]
Gamma   = 1.5;              % FET channel noise factor
g_m     = 0.030;            % FET transconductance [Siemens] [mS]
I_2     = 0.562;            % Noise BW factor
I_3     = 0.0868;           % Noise BW factor

% Room size
x_max       = 5;            % Room size x-axis         % [SCA+2022]
y_max       = 5;            % Room size y-axis         % [SCA+2022]
z_max       = 3;            % Room size z-axis         % 3 in [SCA+2022]

m = -(log(2)/log(cosd(Psi)));
G = (a^2)/((sind(Phi_FoV)^2));

save("config.mat")
end