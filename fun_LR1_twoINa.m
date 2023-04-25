function [G, Iion, Ivec, dX, Iall] = fun_LR1_twoINa(t,X,p,S)
if ~isfield(p,'dt')
    p.dt = 1;
end
if ~isfield(p,'vclamp_flag')
    p.vclamp_flag = 0;
end
dt = p.dt;
Nm = p.Npatches;

Vm = X(1:Nm);
if p.vclamp_flag
    Vm = get_voltage_clamp(t, p.t_clamp, p.V_clamp); 
end

G = X(Nm+1:end);
dX = zeros(length(Vm)+length(G),1);

M_1 = X(Nm+1:2*Nm); % order has to be consistent with the ICs
H_1 = X(2*Nm+1:3*Nm);
J_1 = X(3*Nm+1:4*Nm);
D = X(4*Nm+1:5*Nm);
F = X(5*Nm+1:6*Nm);
Xi = X(6*Nm+1:7*Nm);
Ca_i = X(7*Nm+1:8*Nm);
M_2 = X(8*Nm+1:9*Nm); % order has to be consistent with the ICs
H_2 = X(9*Nm+1:10*Nm);
J_2 = X(10*Nm+1:11*Nm);

Nao = S(1:Nm);   % extracellular Na
Ko = S(Nm+1:2*Nm); % extracellular Ko
Cao = S(2*Nm+1:3*Nm); % extracellular Ca


[INa_1,am_1,bm_1,ah_1,bh_1,aj_1,bj_1] = comp_INa1(Vm,M_1,H_1,J_1,Nao, p);
[INa_2,am_2,bm_2,ah_2,bh_2,aj_2,bj_2] = comp_INa2(Vm,M_2,H_2,J_2,Nao, p);
[Isi,ad,bd,af,bf] = comp_Isi91(Vm,D,F,Ca_i,Cao,p);
[IK,ax,bx] = comp_IK(Vm,Xi,Nao,Ko,p);
IK1 = comp_IK1(Vm,Ko,p);
IKp = comp_IKp(Vm,Ko,p);
Ib = p.f_I(p.iib,:)'.*p.Ib_max.*(Vm+59.87);

Ivec = nan(Nm*3,1);    % ionic currents
Ivec(1:Nm) = INa_1+INa_2;      % Na+
Ivec(Nm+1:2*Nm) = IK + IK1 + IKp; % K+
Ivec(2*Nm+1:3*Nm) = Isi;      % Ca2+


taum_1 = 1./(am_1+bm_1); minf_1 = am_1.*taum_1;
tauh_1 = 1./(ah_1+bh_1); hinf_1 = ah_1.*tauh_1;
tauj_1 = 1./(aj_1+bj_1); jinf_1 = aj_1.*tauj_1;
taum_2 = 1./(am_2+bm_2); minf_2 = am_2.*taum_2;
tauh_2 = 1./(ah_2+bh_2); hinf_2 = ah_2.*tauh_2;
tauj_2 = 1./(aj_2+bj_2); jinf_2 = aj_2.*tauj_2;
taud = 1./(ad+bd); dinf = ad.*taud;
tauf = 1./(af+bf); finf = af.*tauf;
taux = 1./(ax+bx); xinf = ax.*taux;

if p.mLR1   % modified LR1, as in Fenton-Karma 1998
    taud = taud/2;
    tauf = tauf/2;
end

% Rush-Larson for HH gates
G(1:Nm) = minf_1 - (minf_1 - M_1).*exp(-dt./taum_1);       % M
G(Nm+1:2*Nm) = hinf_1 - (hinf_1 - H_1).*exp(-dt./tauh_1);         % H
G(2*Nm+1:3*Nm) = jinf_1 - (jinf_1 - J_1).*exp(-dt./tauj_1);       % J
G(3*Nm+1:4*Nm) = dinf - (dinf - D).*exp(-dt./taud);       % D
G(4*Nm+1:5*Nm) = finf - (finf - F).*exp(-dt./tauf);       % F
G(5*Nm+1:6*Nm) = xinf - (xinf - Xi).*exp(-dt./taux);  
G(6*Nm+1:7*Nm) = Ca_i + dt*(-1e-4*Isi + 0.07*(1e-4-Ca_i)); % Ca_i
G(7*Nm+1:8*Nm) = minf_2 - (minf_2 - M_2).*exp(-dt./taum_2);       % M
G(8*Nm+1:9*Nm) = hinf_2 - (hinf_2 - H_2).*exp(-dt./tauh_2);         % H
G(9*Nm+1:10*Nm) = jinf_2 - (jinf_2 - J_2).*exp(-dt./tauj_2);       % J


Istim = p.stim_amp*(mod(t,p.bcl)<p.stim_dur).*p.indstim;
Iion = -(Istim - INa_1 - INa_2 - Isi - IK - IK1 - IKp - Ib);  % uA

% all ionic currents
[Ncurrents,~] = size(p.f_I);
Iall = nan(Nm*Ncurrents,1);

Iall(p.iina1:Ncurrents:end) = INa_1;
Iall(p.iina2:Ncurrents:end) = INa_2;
Iall(p.iisi:Ncurrents:end) = Isi;
Iall(p.iik:Ncurrents:end) = IK;
Iall(p.iik1:Ncurrents:end) = IK1;
Iall(p.iikp:Ncurrents:end) = IKp;
Iall(p.iib:Ncurrents:end) = Ib;

if Nm == 1
    dX(1:Nm) = -Iion/p.Ctot;
    dX(Nm+1:2*Nm) = (minf_1 - M_1)./taum_1;       % M
    dX(2*Nm+1:3*Nm) = (hinf_1 - H_1)./tauh_1;         % H
    dX(3*Nm+1:4*Nm) = (jinf_1 - J_1)./tauj_1;       % J
    dX(4*Nm+1:5*Nm) = (dinf - D)./taud;       % D
    dX(5*Nm+1:6*Nm) = (finf - F)./tauf;       % F
    dX(6*Nm+1:7*Nm) = (xinf - Xi)./taux;
    dX(7*Nm+1:8*Nm) = (-1e-4*Isi + 0.07*(1e-4-Ca_i)); % Ca_i
    dX(8*Nm+1:9*Nm) = (minf_2 - M_2)./taum_2;       % M
    dX(9*Nm+1:10*Nm) = (hinf_2 - H_2)./tauh_2;         % H
    dX(10*Nm+1:11*Nm) = (jinf_2 - J_2)./tauj_2;       % J
end
end

function [INa,alpham,betam,alphah,betah,alphaj,betaj] = comp_INa1(V,m,H,J,Nao,data)
% INa    the fast sodium current in mammalian ventricular cells
%

ENa =data.RTF*log(Nao./data.Na_i);       % Nernst potential of Na, mV
INa = data.INa_max.*data.f_I(data.iina1,:)'.*m.*m.*m.*H.*J.*(V-ENa);

alpham = am(V); betam = bm(V);
alphah = ah(V); betah = bh(V);
alphaj = aj(V); betaj = bj(V);
end

function [INa,alpham,betam,alphah,betah,alphaj,betaj] = comp_INa2(V,m,H,J,Nao,data)
% INa    the fast sodium current in mammalian ventricular cells
%

ENa =data.RTF*log(Nao./data.Na_i);       % Nernst potential of Na, mV
INa = data.INa_max.*data.f_I(data.iina2,:)'.*m.*m.*m.*H.*J.*(V-ENa);

alpham = am(V+data.Vact_shift); betam = bm(V+data.Vact_shift);
alphah = ah(V+data.Vinact_shift); betah = bh(V+data.Vinact_shift);
alphaj = aj(V+data.Vinact_shift); betaj = bj(V+data.Vinact_shift);
end


function ah = ah(V)
a=1-1./(1+exp(-(V+40)/0.24));
ah= a.*0.135.*exp((80+V)./(-6.8));
%  ah= a.*0.135.*exp((80+V)./(-15.8));

end

function bh = bh(V)
a=1-1./(1+exp(-(V+40)/0.24));
bh= (1-a)./(0.13*(1+exp((V+10.66)/(-11.1)))) +(a).*(3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V));
% bh= (1-a)./(0.13*(1+exp((V+10.66)/(-11.1)))) +(a).*(.1*3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V));

end

function aj = aj(V)
a=1-1./(1+exp(-(V+40)/0.24));
aj =  a.*(-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V)).*(V+37.78)./(1+exp(0.311*(V+79.23)));
end

function bj = bj(V)
a=1-1./(1+exp(-(V+40)/0.24));
bj= (1-a).*(0.3*exp(-2.535e-7*V)./(1+exp(-0.1*(V+32))))+(a).*(0.1212*exp(-0.01052*V)./(1+exp(-0.1378*(V+40.14))));
end

function am = am(V)
am = 0.32*(V+47.13)./(1-exp(-0.1*(V+47.13)));
end

function bm = bm(V)
bm = 0.08*exp(-V/11);
end

function [IK,alphax,betax] = comp_IK(V,x,Nao, Ko,data)
% IK    Time-dependent potassium current
Xi = xi(V);
EK = data.RTF.*log((Ko+data.PNa_K*Nao)/(data.K_i+data.PNa_K*data.Na_i));        % mV
sq=sqrt(Ko/5.4); 
GK = data.IK_max.*sq;         % mS/um^2

IK = data.f_I(data.iik,:)'.*GK.*Xi.*x.*(V-EK);
alphax = ax(V); betax = bx(V);
end

function Xi = xi(V)
Xi = (V<=-100)+ ...
    (V > -100).*(2.837*(exp(0.04*(V+77))-1)./((V+77).*exp(0.04*(V+35))));
end

function ax = ax(V)
ax = 0.0005*exp(0.083*(V+50))./(1+exp(0.057*(V+50)));
end

function bx = bx(V)
bx = 0.0013*exp(-0.06*(V+20))./(1+exp(-0.04*(V+20)));
end 

function IK1 = comp_IK1(V,Ko,data)
% IK1    Time-independent potassium current

EK1 = data.RTF*log(Ko./data.K_i);
ak1 = 1.02./(1+exp(0.2385*(V-EK1-59.215)));
bk1 = (0.49124*exp(0.08032*(V-EK1+5.476))+exp(0.06175*(V-EK1-594.31)))./...
    (1+exp(-0.5143*(V-EK1+4.753)));
sq=sqrt(Ko/5.4); 

gK1 =data.IK1_max.*sq.*ak1./(ak1+bk1);
IK1 = data.f_I(data.iik1,:)'.*gK1.*(V-EK1);
end

function [Is,alphad,betad,alphaf,betaf] = comp_Isi91(V,d,f,Ca_i,Cao,data)
% Isi    the slow inward current
% Esi = 7.7 - 13.0287*log(Ca_i);
Esi = data.RTF/2*log(Cao./Ca_i);
Is =  data.f_I(data.iisi,:)'.*data.Isi_max.*d.*f.*(V-Esi);

alphad = ad(V); betad = bd(V);
alphaf = af(V); betaf = bf(V);
end

function ad = ad(V)
ad = 0.095*exp(-0.01*(V-5))./(1+exp(-0.072*(V-5)));
end

function bd = bd(V)
bd = 0.07*exp(-0.017*(V+44))./(1+exp(0.05*(V+44)));
end

function af = af(V)
af = 0.012*exp(-0.008*(V+28))./(1+exp(0.15*(V+28)));
end

function bf = bf(V)
bf = 0.0065*exp(-0.02*(V+30))./(1+exp(-0.2*(V+30)));
end

function IKp = comp_IKp(V,Ko,data)
EK1 = data.RTF*log(Ko./data.K_i);
IKp = data.f_I(data.iikp,:)'.*data.IKp_max.*(V-EK1)./(1+exp((7.488-V)/5.98));
end


