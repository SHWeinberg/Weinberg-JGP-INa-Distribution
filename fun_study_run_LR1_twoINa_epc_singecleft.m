function [cv, tup, ts, phi_i, phi_cleft, phi_L, phi_R, INa1, INa2, Na_cleft] = ...
    fun_study_run_LR1_twoINa_epc_singecleft(Ncell, D, fNa1, fNa2, ...
    w, locINa1, locINa2, L, r, total_fNa, single_beat_flag, bcl, nbeats)

% inputs:
% Ncell: # of cells
% D: diffusion coefficient (cm^2/s), used to calculate ggap (along with r, L)
% fNa1, fNa2: scaling factors for two INa subpopulations (should sum to 1)
% w: cleft width, um
% locINa1, locINa2: ID localization for two INa subpoplations (between 0 and 1)
% L: cell length, um
% r: cell radius, um
% total_fNa: scaling factor for total INa conductance
% single_beat_flag: flag to stop sim after wave propagated to end of cable

% time parameters
% bcl = 500;      % basic cycle length ms
% nbeats = 1;      % number of beats
T = bcl*nbeats;

single_cell_ic_flag = 1; % flag to run single cell sim to determine cable initial conditions 
if nargout > 1
    out_val = 10; % use smaller sampling time step, if outputing state variables
else
    out_val = 1;
end

% estimate propagation time down cable, based on a set conduction velocity value,
% used to determine twin (see below)
cv_est = 10*1e3/100;  % um/ms (D = 0.1 ~ 20 cm/s, use 10 cm/s)
Lcable = L*Ncell; % um
tprop_est = Lcable/cv_est;  % ms

% model-independent parameters
if total_fNa < 5
    dt_scale = 1;
else % smaller time step needed for large total_fNa
    dt_scale = 0.1;
end

% different numerical integration and sampling methods
if 0
    dtparams.method = 'fixed';  % 'fixed', 'adapt', or 'adapt_samp'
    dtparams.dt1 = dt_scale*1e-2;      % ms
    dtparams.dt1_samp = dtparams.dt1*1;
    dtparams.dt2 = dtparams.dt1*50;
    dtparams.dt2_samp = dtparams.dt2*1;
    dtparams.dtS1 = dtparams.dt1*0.2; % time step for cleft concentration integration
    dtparams.dtS2 = dtparams.dt1*1;
    dtparams.twin = max(25, round(tprop_est)+10); % defines window using dt1
elseif 0
    dtparams.method = 'adapt';
    dtparams.dt_min = dt_scale*1e-2;  % min time step used during upstroke (< twin), ms
    dtparams.dt_max = dtparams.dt_min*50;  % max time step, ms
    dtparams.dtS_min = dtparams.dt_min*0.2; % min time step for cleft concentration integration, used while dt = dt_min, should be <= dt_min
    dtparams.dtS_max = dtparams.dt_min*2;  % max time step for cleft concentration, must be <= dt_max
    dtparams.dV_max = 0.5;  % max voltage change, mV
    dtparams.dS_max = 0.01;  % max cleft concentration change, mM
    dtparams.twin = max(25, round(tprop_est)+10); % defines window using dt_min
else
    dtparams.method = 'adapt_samp';
    dtparams.dt_min = dt_scale*1e-2;  % min time step used during upstroke (< twin), ms
    dtparams.dt_min_samp = dtparams.dt_min*10/out_val; % sampling
    dtparams.dt_max = dtparams.dt_min*1;  % max time step, ms
    dtparams.dt_max_samp = dtparams.dt_max*1;
    dtparams.dtS_min = dtparams.dt_min*0.2; % min time step for cleft concentration integration, used while dt = dt_min, should be <= dt_min
    dtparams.dtS_max = dtparams.dt_min*2;  % max time step for cleft concentration, must be <= dt_max
    dtparams.dV_max = 0.5;  % max voltage change, mV
    dtparams.dS_max = 0.01;  % max cleft concentration change, mM
    dtparams.twin = max(25, round(tprop_est)+10); % defines window using dt_min
end
Vthresh = -60;  % activation/repol threshold, mV
% single_beat_flag = 1; % flag to stop sim after single beat propagated to end of cable
trange = [0 T];
% trange = [max(0,(nbeats-4)*bcl-10) nbeats*bcl];  % time range to output values
% trange = [0 nbeats*bcl];  % time range to output values

Mdisc = length(w);      % number of patches at ID
p_ext = 150*10;  % extracellular resistivity, k-ohm*um

f_bulk = 1;
Gb_mat = 8*pi*f_bulk*w/p_ext; % cleft-bulk conductance, mS
Gc_array = 0;

% cell geometry
Aax = 2*pi*r*L; % patch surface area, um^2
IDarea_vec = pi*r^2/Mdisc*ones(Mdisc,1);  % individual ID membrane patch surface area, um^2
Ad = sum(IDarea_vec);
Atot = 2*Ad + Aax;  % total surface area, um^2

fVol_cleft = 1; % scaling factor
Vol_cleft_vec = fVol_cleft*IDarea_vec.*w; % cleft volume, um^3

% bulk extracellular concentrations, baseline values
Ko = 5.4;                  % mM, 5.4
Nao = 140;                 % mM, 140
Cao = 1.8;                 % mM, 1.8

% ID localization (vs axial patch)
locUniform = 2*Ad/Atot;  % value for uniform distribution (proportional to area)

Npatches = 2*Mdisc*(Ncell-1) + Ncell; % total number of membrane patches

% storage parameters
storeparams.int_store_all = 100;  % ms, interval for storing all state variables
% gating variable indices to store (model specific), [ ] = store all variables, NaN = store no variables
storeparams.indG_store = [1 2 7];

[p, x0] = InitialConstants_LR91_twoINa(Atot);

% order is determined by code in fun_name
% INa1, Isi, IK, IK1, IKp, Ib, INa2
Ncurrents = 7;
scaleI = ones(1,Ncurrents);

p.iina1 = 1; p.iisi = 2; p.iik = 3;
p.iik1 = 4; p.iikp = 5; p.iib = 6;
p.iina2 = 7;

loc_vec = zeros(1, Ncurrents); % ID localization vec
loc_vec(p.iina1) = locINa1;
loc_vec(p.iina2) = locINa2;

scaleI(p.iina1) = fNa1;
scaleI(p.iina2) = fNa2;

% estimated from Lin,..., Delmar, Heart Rhythm 2011
% note the difference in sign (-/+) with the JGP manuscript, due to how the
% equations are defined in the code, i.e., (V-Vshift) vs (V+Vshift)
p.Vact_shift = 5.5; % mV, shift in SS activation curve (+ means left shift)
p.Vinact_shift  = -7.8; % mV, shift in SS inactivation curve, (- means right)

% ionic model-specific parameters
ionic_fun_name = 'fun_LR1_twoINa';

% initial conditions
p.Nstate = 8+3;  % number of state variables, including Vm, per patch
p.mLR1 = 1; % flag for modified LR1 model with Ca2+ speedup

% stimulus parameters
p.stim_dur = 1;   % ms
p.stim_amp = 80e-8*Atot;    % uA
stim_cells = 1:5;
p.istim = (stim_cells-1)*(2*Mdisc+1)+1;  % indices of axial membrane patch
p.indstim = ismember((1:Npatches)',p.istim);

loc_mat = loc_vec;

p.L = L; p.r = r;  % um
% extracellular concentrations
p.K_o = Ko;                  % mM
p.Na_o = Nao;                 % mM
p.Ca_o = Cao;                 % mM
% cleft ionic concentration parameters
p.cleft_conc = [p.Na_o; p.K_o; p.Ca_o];
clamp_flag = [1; 1; 1];  % clamp cleft ionic concentrations Na, K, Ca,
% 0 = not clamped, 1 = clamped (fixed)
zvec = [1; 1; 2];  % charge valence

if single_cell_ic_flag % option to run single cell simulations to obtain new initial conditions, fixed extracellular ionic concentrations
    nbeats_single = 10;
    % use fixed two time step method for single cell sims
    dt1 = 1e-2; dt2 = 1e-1;
    psingle = p; psingle.dt = dt1; psingle.bcl = bcl;
    psingle.Npatches = 1; psingle.f_I = scaleI';
    psingle.indstim = 1;
    psingle.celltype = 0;
    psingle.ind_tau_ip = 1; psingle.ind_tau_im = 1; psingle.tau_Nai_mat = 1; psingle.tau_Cai_mat = 1; psingle.tau_Ki_mat = 1;
    Cm = 1*1e-8;      % membrane capacitance, uF/um^2
    psingle.Ctot = (2*sum(IDarea_vec) + 2*pi*r*L)*Cm;   % total cell capacitance, uF
    Vm_single(1) = x0(1); Gsingle = x0(2:end)';
    ionic_fun = str2func(['@(t,x,p,S) ',ionic_fun_name,'(t,x,p,S)']);

%     Vvec = Vm_single; tvec = 0;
    options = odeset('MaxStep',1,'InitialStep',2e-2);
    tic;
    [~,X] = ode15s(@(t,x) get_ionic_fun_dX(ionic_fun,t,x,psingle, p.cleft_conc),...
        0:dt2:bcl*nbeats_single,[Vm_single;Gsingle], options);

    x0 = X(end,:);

end

p.phii_0 = x0(1);
g0_vec = zeros(Npatches*(p.Nstate-1),1);
for i = 1:p.Nstate-1
    g0_vec(1 + Npatches*(i-1):Npatches*i) = x0(i+1);
end
p.g0_vec = g0_vec;

[ts, phi_mat, Gmat, Smat, Iall_mat, tup, trepol, pout, cv, store_out] =  generic_IDnano_epc_v5(ionic_fun_name, ...
    p, Ncell, Mdisc, Gb_mat, Gc_array, D, bcl, loc_mat, IDarea_vec, scaleI, Vol_cleft_vec, clamp_flag, ...
    zvec, dtparams, T, Vthresh, single_beat_flag, trange, storeparams);

if nargout > 1
    Nphi = 3*Ncell-2 + Mdisc*(Ncell-1);  % number of voltage state variables
    Ncleft_comp = (Ncell-1)*Mdisc;  % number of cleft compartments
    Npatches = 2*Mdisc*(Ncell-1) + Ncell; % number of membrane patches

    Na_cleft = Smat(1:Ncleft_comp,:);

    ind_Vm = ((1:Ncell)-1)*(2*Mdisc+1)+1;  % indices of axial membrane patch
    ind_phim = ((1:Ncell)-1)*(Mdisc+3)+1; % indices of axial potentials
    ind_phiR = ((1:Ncell-1)-1)*(Mdisc+3)+2; % indices of right potentail
    ind_phiL =  ((2:Ncell)-1)*(Mdisc+3);  % indices of left potentials
    ind_phic = setdiff(1:Nphi, [ind_phim ind_phiR ind_phiL]);  % indices of cleft potentials
    ind_VR = []; ind_VL = [];
    for i = 1:Ncell-1, for j = 1:Mdisc, ind_VR = [ind_VR, (i-1)*(2*Mdisc+1)+j+1]; end, end
    for i = 2:Ncell, for j = 1:Mdisc, ind_VL = [ind_VL, (i-2)*(2*Mdisc+1)+Mdisc+j+1]; end, end

    icell = round(Ncell/2);  % indices for a specific cell
    i_phii_cell = icell;
    i_VL_cell = (icell-2)*Mdisc+1:(icell-2)*Mdisc+Mdisc;
    i_VR_cell = (icell-1)*Mdisc+1:(icell-1)*Mdisc+Mdisc;

    Vm_mat = constructVm_mat(phi_mat, Npatches, Mdisc, Ncell);

    phi_i = phi_mat(ind_phim,:);
    phi_cleft = phi_mat(ind_phic,:);
    phi_L = phi_mat(ind_phiL,:);
    phi_R = phi_mat(ind_phiR,:);

    % note: Vm = phi_i = Vm_mat(ind_Vm,:)
    VL = Vm_mat(ind_VL,:);
    VR = Vm_mat(ind_VR,:);

    Ncurrents = length(scaleI);
    INa1 = Iall_mat(pout.iina1:Ncurrents:end,:);
    INa2 = Iall_mat(pout.iina2:Ncurrents:end,:);

end