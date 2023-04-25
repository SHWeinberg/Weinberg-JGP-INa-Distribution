function ggap = convert_D_to_ggap(D, r, L)


% convert D (cm^2/s) to ggap (mS) based on cell dimensions, capacitance

Aax = 2*pi*r*L; % patch surface area, um^2
Ad = pi*r^2;    % disc surface area, um^2
Atot = 2*Ad + Aax;  % total surface area, um^2
Cm = 1*1e-8;      % membrane capacitance, uF/um^2
Ctot = Atot*Cm; % total capacitance, uF

dx = L/1000;    % mm
Dnorm = D*100/1000;  % mm^2/ms
R = dx^2./Dnorm;  % normalized resistance (ms)
Rgap = R./Ctot;  % k-ohms ( = ms/uF)
ggap = 1./Rgap;  % mS, gap junctional conductance
