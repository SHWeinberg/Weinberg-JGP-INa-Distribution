function [data, x0] = InitialConstants_LR91_twoINa(Atot)

% Atot - total cell area in um^2

% Standard ionic concentrations for ventricular cells
data.K_o = 5.4;                  % mM
data.K_i = 145;                  % mM
data.Na_o = 140;                 % mM
data.Ca_o = 1.8;                 % mM
data.Na_i = 10;                  % mM
data.F = 96.5;                   % Faraday constant, coulombs/mmol
data.R = 8.314;                  % gas constant, J/K
data.T = 273+37;                 % absolute temperature, K 
data.RTF=(data.R*data.T/data.F); % mV
% data.sqrt=sqrt(data.K_o/5.4); 

data.PK = 1.66e-6;                 % permability of K 
data.PNa_K = 0.01833;              % permability ratio of Na to K

% factor of 1e-8 converts from mS/cm^2 to mS/um^2 to mS
data.Isi_max=0.09e-8*Atot;               % mS
data.IKp_max=0.0183e-8*Atot;
data.Ib_max=0.03921e-8*Atot;
data.INa_max=16e-8*Atot;
data.IK1_max=0.6047e-8*Atot;
data.IK_max=0.282e-8*Atot;


%the initial conditions
% % Initial Gate Conditions */

v_init=-84.5286 ;  % mV
m_init = 0.0017; % sodium current activation gate
h_init =   0.9832;  % sodium current fast inactivation
J_init = 0.995484;   %  slow inactivation
d_init = 0.000003; % Calcium activation gate
f_init =  1.0000 ;  % Calcium  inactivation gate
X_init =  0.0057 ; % activation gate
Ca_init = 0.0001; % Calcium, mM


x0=[ v_init m_init h_init J_init d_init f_init X_init Ca_init m_init h_init J_init];
end
