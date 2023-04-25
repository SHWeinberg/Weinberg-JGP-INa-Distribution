function [ts, phi_mat, Gmat, Smat, Iall_mat, tup, trepol, p, cv_est, store_out] =  ...
    generic_IDnano_epc_v5(ionic_fun_name, p, Ncell, Mdisc, Gb_mat, Gc_array,...
    D, bcl, loc_mat, IDarea_vec, scaleI, Vol_cleft_vec, clamp_flag, zvec,...
    dtparams, T, Vthresh, single_beat_flag, trange, storeparams)
%
% outputs:
% 
% ts: time vector (ms), Ns x 1 vector, values at specified sampling times (set in dtparams) 
%
% phi_mat: potential matrix (mV), Nphi x Ns matrix
%
% Gmat: gating variable matrix, Npatches*(Nstate-1) x Ns matrix (units
% as defined in ionic model)
%
% Smat: cleft ion concontration matrix (mM), 3*Ncleft_comp x Ns matrix
%
% Iall_mat: all ionic currents matrix (uA), Ncurrents*Npatches x Ns (order
% as defined in ionic model): grouped by membrane patches, 
% [Current 1 (Patch 1),..., Current Nc (Patch 1),
%  Current 1 (Patch 2),..., Current Nc (Patch 2),..., 
%  Current 1 (Patch Np),..., Current Nc (Patch Np)]
% 
% tup, trepol: matrix of upstroke and repolarization times (ms), 
% Npatches x nbeats (typically - could differ if conduction block, extra beats, etc)
%
% p: output structure of defined parameters
%
% cv_est: nbeats x 1 vector of estimated conduction velocity values (cm/s)
%
% store_out: structure consisting of all values for phi, G, S at specified
% time intervals (set in storeparams)
%
% inputs:
%
% ionic_fun_name: function name for ionic model, specific input order (see
% fun_LR1 for example)
%
% p: structure containing model-specific parameters
%
% Ncell: # of cells
%
% Mdisc: # of ID membrane patches/cleft compartments
%
% Gb_mat: bulk conductance matrix (mS), Ncell x Mdisc matrix (or 1 x Mdisc
% if assuming all cells are the same),
% such that for cleft i (right of cell i), Gb(i,j) = conductance between patch j and bulk
%
% Gc_array: cleft conductance array (mS), Mdisc x Mdisc x Ncell array (or Mdisc x Mdisc matrix
% if assuming all cells are the same)
% such that Gc(j,k,i) = conductance between patch j and k for cleft i
% Note: first 2 dimensions should be symmetric, i.e., Gc(j,k,i) = G(k,j,i) 
% but must be specified as such
%
% D: diffusion coefficient (cm^2/s), Ncell-1 x 1 vector (or scalar if all gap
% junctions the same), used to calculate gap junction conductance
%
% bcl: basic cycle length (ms)
%
% loc_mat: ID localization matrix [0,1], 
% if loc_mat Mdisc x Ncurrents (where Ncurrents = # of currents in ionic model)
% such that loc_mat(j,q) = fraction of total cell conductance for current q
% in ID patch j.  Assumed to be the same for all IDs.
% Note: loc_ID = sum(loc_mat), such that loc_ID(q) = fraction of total cell
% conductance for current q at the ID, and loc_ax = 1-loc_ID, such that
% loc_ax(q) = fraction of total cell conductance for current q on the axial
% membrane.  Thus, sum(loc_mat) + loc_ax = 1
% 
% if loc_mat 2*Mdisc x Ncurrents, then first Mdisc rows are the
% pre-junctional localization, and the second Mdisc rows are the
% post-junctional localization
%
% IDarea_vec: ID membrane patch area (um^2), Mdisc x 1.  Assumed to be
% symmetric (both sides of cleft) and same for all cells
%
% scaleI: ionic current scaling factor vector, 1 x Ncurrents vector
% or 1 x Ncurrents cell array, where scaleI{q} can either be a scalar or 
% Ncell x 1 vector specifying scaling factor for each cell
%
% Vol_cleft_vec: cleft compartment volume vector (um^3), Mdisc x 1
%
% clamp_flag: flag to clamp cleft concentrations (1 = clamped), 1 x 3 (for
% Na, K, Ca)
%
% zvec: charge valence, 1 x 3 (for Na, K, Ca)
%
% dtparams structure:
% method: 'fixed' or 'adapt' for fixed or adaptive time steps
% Params for 'fixed' method
% dt1: time step (msec) - used from time of stimulus to twin (time period 1)
% dt1_samp: sampling time interval (for time period 1)
% dt2: time step (msec) - used from twin to time of next stimulus (time period 2)
% dt2_samp: sampling time interval (for time period 2)
% Note: by design dti_samp >= dti for both i = 1, 2
% dtS1: time step (msec) - used for cleft concentration integration (time period 1)
% dtS2: time step (msec) - used for cleft concentration integration (time period 2)
% twin: duration to utilize time step dt1.  Use values such that dt1<=dt2
%
% Params for 'adapt' method
% dt_min: min time step (msec) - minimum time step, always used from time of stimulus to twin 
% dt_max: max time step (msec) - maximum time step, used from twin to time
% of next stimulus
% dV_max: max voltage change (mV) - maximum value for voltage change (all
% potential values) 
% dtS_min: min time step (msec) - used for cleft concentration integration
% when dt = dt_min
% dtS_max: max time step (msec) - used for cleft concentration integration
% dS_max: max cleft concentration change (mM) - maximum value for cleft
% concentration change (all values)
% twin: duration to utilize time step dt_min.  
%
% Params for 'adapt_samp' method: 
% Same as 'adapt' with addition parameters for sampling
% dt_min_samp: sampling time interval between time of stimulus and twin
% dt_max_samp: sampling time interval between twin and time of next stimulus
%
% T: total simulation time (msec)
% Vthresh: voltage threshold for activation/repolarization (mV)
% single_beat_flag: flag to only simulate first activation wave
% trange: time interval to save state variables (msec), 1 x 2 vector
%
% storeparams structure:
% int_store_all: time interval (ms) for storing all state variables
% indG_store: indices for variables in Gmat (model-dependent) to store
% every sampling time, 
% empty vector [ ] = store all variables, NaN = store no variables
%
% Units:
%
% voltage in mV
% current in uA
% conductance in mS (= 1/k-ohm)
% resistance in k-ohm (= 1/mS)
% capacitance in uF    (Recall: uF*mV/ms = uA)
% time in ms
% length in um
% concentration in mM
%
% New in v5:
% 1. localization can differ between pre- and post-junctional membranes. If
% loc_mat has 2*Mdisc rows, then the first Mdisc rows define pre-junctional
% localization, and the second Mdisc rows define post-junctional
% localization
%
% New in v4:
% 1. storeparams structure enables storing specific gating variables in
% output matrix Gmat, and stores all state variables at specified intervals
% 2. output Iall_array, all ionic currents
% 3. scaleI can be a cell array to specify spatially varying ion
% conductance scaling factors
% 
%
% New in v3:
% 1. parameter structure p is an output (used to calculate
% currents, etc for analysis)
% 2. use of adaptive time steps, one for immediately following stimulus, one for
% the duration of beat, or an adaptive time step method

if nargin == 0
    Ncell = 20;     % number of cells
    bcl = 400;      % basic cycle length ms
    nbeats = 2;
    % estimate Gb and Gc from cleft widths
     %   w_vec = [10e-3];      % Mdisc x 1, cleft widths vector, um
     w_vec = 1e-3*randi([10,50],4,1);
    Mdisc = length(w_vec);      % number of patches at ID
    p_ext = 150*10;  % extracellular resistivity, k-ohm*um
    
    % Rradial = p_ext/(8*pi*w);   % radial cleft resistance, k-ohm
        %     Gb_mat =   8*pi*repmat(w_vec',Ncell,1)/p_ext; % mS, bulk conductance, Ncell x Mdisc matrix
    f_disc = 1;
    f_bulk = 1;
        Gb_mat = 8*pi*f_bulk*w_vec'/p_ext; % mS, bulk conductance, 1 x Mdisc vector
        % for cleft i (right of cell i), Gb(i,j) = conductance between patch j and bulk

        Gci = zeros(Mdisc, Mdisc);
        for j = 1:Mdisc
            for k = j+1:Mdisc
                if j ~= k      % random connectivity
                    Gci(j,k) = (rand<0.5)*f_disc*8*pi*(w_vec(j) + w_vec(k))/2/p_ext;  % mS
                    Gci(k,j) = Gci(j,k);
                end
            end
        end
        Gc_array = Gci;  % Mdisc x Mdisc matrix (assumes all cells are identical)
        %     Gc_array = nan(Mdisc,Mdisc,Ncell);
        % mS, cleft conductances, M x M x N array, Gc(j,k,i) = conductance between patch j and k for cleft i
        %     for i = 1:Ncell, Gc_array(:,:,i) = Gci; end
    
    
    loc_vec = [0.9 0 0 0.5 0 0];  % ID localization vec
    % order is determined by code in fun_name
    % INa, Isi, IK, IK1, IKp, Ib
    scaleI = [1 1 1 1 1 1];

    
    % equal distribution in the ID nanodomains
    loc_mat = repmat(loc_vec, Mdisc, 1)/Mdisc;  % Mdisc x Ncurrents, ID localization matrix
    
    
    D = 1*ones(Ncell-1,1);      % diffusion coefficient, cm^2/s, used to determine ggap
    D(round(Ncell/2)) = 1;
    
    dtparams.method = 'adapt_samp';
    dtparams.dt_min = 1e-2;  % min time step used during upstroke (< twin), ms
    dtparams.dt_min_samp = dtparams.dt_min*1; % sam
    dtparams.dt_max = dtparams.dt_min*10;  % max time step, ms
    dtparams.dt_max_samp = dtparams.dt_max*10;
    dtparams.dtS_min = dtparams.dt_min*0.2; % min time step for cleft concentration integration, used while dt = dt_min, should be <= dt_min
    dtparams.dtS_max = dtparams.dt_min*2;  % max time step for cleft concentration, must be <= dt_max
    dtparams.dV_max = 0.5;  % max voltage change, mV
    dtparams.dS_max = 0.01;  % max cleft concentration change, mM
    dtparams.twin = 25; % defines window using dt_min
    % dtparams.method = 'fixed';
%     dtparams.dt1 = 1e-2;      % ms
%     dtparams.dt1_samp = dtparams.dt1;
%     dtparams.dt2 = 5e-1;
%     dtparams.dt2_samp = dtparams.dt2;
%     dtparams.dtS1 = dtparams.dt1*0.2;
%     dtparams.dtS2 = dtparams.dt1;
%     dtparams.twin = 25;  % time interval for time step of dt1
    T = bcl*nbeats;
    Vthresh = -60;  % activation/repol threshold, mV
    single_beat_flag = 1;
    
    % cell geometry parameters
    L = 100;        % cell length, um
    r = 11;         % cell radius, um
    Aax = 2*pi*r*L; % patch surface area, um^2
    %     Ad = pi*r^2;    % disc surface area, um^2
    IDarea_vec = pi*r^2/Mdisc*ones(Mdisc,1);  % individual ID membrane patch surface area, um^2
    Ad = sum(IDarea_vec);
    Atot = 2*Ad + Aax;  % total surface area, um^2
    
    fVol_cleft = 1;
    Vol_cleft_vec = fVol_cleft*IDarea_vec.*w_vec; % cleft volume, um^3
    
    
    % ionic model-specific parameters
    ionic_fun_name = 'fun_LR1';
    
    [p, x0] = InitialConstants_LR91(Atot);
    p.iina = 1; p.iisi = 2; p.iik = 3;
    p.iik1 = 4; p.iikp = 5; p.iib = 6;
    p.Nstate = 8;
    p.phii_0 = x0(1);
    Npatches = 2*Mdisc*(Ncell-1) + Ncell; % number of membrane patches
    
    g0_vec = zeros(Npatches*(p.Nstate-1),1);
    for i = 1:p.Nstate-1
        g0_vec(1 + Npatches*(i-1):Npatches*i) = x0(i+1);
    end
    p.g0_vec = g0_vec;
    p.mLR1 = 1; % flag for modified LR1 model with Ca2+ speedup
    
    % stimulus parameters
    p.stim_dur = 1;   % ms
    p.stim_amp = 80e-8*Atot;    % uA
    p.istim = ((1:5)-1)*(2*Mdisc+1)+1;  % indices of axial membrane patch
    p.indstim = ismember((1:Npatches)',p.istim);
    
    p.L = L; p.r = r;
    
    % cleft ionic concentration parameters
    p.cleft_conc = [p.Na_o; p.K_o; p.Ca_o];
    clamp_flag = [0; 1; 1];  % Na, K, Ca
    zvec = [1; 1; 2];  % charge valence
    
    trange = [max(0,(nbeats-2)*bcl-10) nbeats*bcl];
    
    storeparams.int_store_all = 100;  % ms
    storeparams.indG_store = [5]; % gating variable indices to store, [ ] = store all variables, NaN = store no variables
    
end

if isempty(trange)
    trange = [0 T];
end
if size(Gb_mat,1)==1
    Gb_mat = repmat(Gb_mat,Ncell,1);
end
if size(Gc_array,3)==1
    Gc_array = repmat(Gc_array,1,1,Ncell);
end
if length(D)==1
    D = D*ones(Ncell-1,1);
end

[tmp_a, ~] = size(loc_mat);
if tmp_a == Mdisc
    loc_matPre = loc_mat/2;
    loc_matPost = loc_mat/2;
else
    loc_matPre = loc_mat(1:Mdisc,:);
    loc_matPost = loc_mat(Mdisc+1:end,:);
end

% t = 0:dt:T;
switch dtparams.method
    case 'fixed'
        dt1 = dtparams.dt1; dt2 = dtparams.dt2;
        dt1_samp = dtparams.dt1_samp; dt2_samp = dtparams.dt2_samp;
        dtS1 = dtparams.dtS1; dtS2 = dtparams.dtS2;
        twin = dtparams.twin;
        ts = get_time_variable(trange, dt1_samp, dt2_samp, twin, bcl);
        Ns = length(ts);
    case 'adapt'
        dt_min = dtparams.dt_min; dt_max = dtparams.dt_max;
        dV_max = dtparams.dV_max;
        dtS_min = dtparams.dtS_min;  dtS_max = dtparams.dtS_max;
        dS_max = dtparams.dS_max;
        twin = dtparams.twin;
        ts_est = get_time_variable(trange, dt_min, dt_max, twin, bcl); % used to estimate length of total state spaces
        Ns = round(1.2*length(ts_est));
    case 'adapt_samp'
        dt_min = dtparams.dt_min; dt_max = dtparams.dt_max;
        dt_min_samp = dtparams.dt_min_samp;
        dt_max_samp = dtparams.dt_max_samp;
        dV_max = dtparams.dV_max;
        dtS_min = dtparams.dtS_min;  dtS_max = dtparams.dtS_max;
        dS_max = dtparams.dS_max;
        twin = dtparams.twin;
        ts = get_time_variable(trange, dt_min_samp, dt_max_samp, twin, bcl);
        Ns = length(ts);     
end
p.bcl = bcl;

L = p.L; r = p.r;   % cell dimensions, um

p_myo = 150*10; % myoplasmic resistivity, k-ohm*um
% p_ext = 150*10;  % extracellular resistivity, k-ohm*um
% Rradial = p_ext/(8*pi*w);   % radial cleft resistance, k-ohm

Cm = 1*1e-8;      % membrane capacitance, uF/um^2
Aax = 2*pi*r*L; % patch surface area, um^2
% Ad = pi*r^2;    % disc surface area, um^2
Ad = sum(IDarea_vec); % disc surface area, um^2
Atot = 2*Ad + Aax;  % total surface area, um^2
Cax = Cm*Aax;    % axial capacitance, uF
Cdisc = Cm*IDarea_vec;  % cleft capacitance, uF, Mdisc x 1 vector

CL = repmat(Cdisc',Ncell, 1); % uF, left side ID patch capacitance, Ncell x Mdisc matrix
CR = repmat(Cdisc',Ncell, 1); % uF, right side ID patch capacitance, Ncell x Mdisc matrix

p.Ctot = Atot*Cm;   % total cell capacitance, uF
p.CL = CL; p.CR = CR; p.Cax = Cax;

Nphi = 3*Ncell-2 + Mdisc*(Ncell-1);  % number of voltage state variables
Npatches = 2*Mdisc*(Ncell-1) + Ncell; % number of membrane patches
Ncleft_comp = (Ncell-1)*Mdisc;  % number of cleft compartments

Vol_cleft_all = repmat(Vol_cleft_vec,  Ncell-1, 1);

p.N = Ncell;
p.M = Mdisc;
p.Npatches = Npatches;
p.Nphi = Nphi;
p.Ncleft_comp = Ncleft_comp;

Rmyo = p_myo*(L/2)/(pi*r^2);    % myoplasmic resistance, k-ohm
gmyo = 1/Rmyo;                  % mS, myoplasmic conductance, scalar

if 0
    fgap = 1; Rgap = 395/fgap;      % gap junction resistance, k-ohm
else
    dx = L/1000;    % mm
    D = D*100/1000;  % mm^2/ms
    R = dx^2./D;  % normalized resistance (ms)
    Rgap = R/p.Ctot;  % k-ohms ( = ms/uF)
end
ggap = 1./Rgap;    % mS, gap junctional conductance, Ncell-1 x 1 vector, between cell i and i+1, for i = 1,...,N-1

switch dtparams.method
    case 'fixed'
        [P1, Q1, C1] = make_IDnano_matrices_v2(Ncell, Mdisc, dt1, gmyo, Cax, CL, CR, ggap, Gb_mat, Gc_array, length(zvec));
        [P2, Q2, C2] = make_IDnano_matrices_v2(Ncell, Mdisc, dt2, gmyo, Cax, CL, CR, ggap, Gb_mat, Gc_array, length(zvec));
        Pinv1 = inv(P1);
        Pinv2 = inv(P2);
    case {'adapt','adapt_samp'}
        [Pmin, Qmin, Cmin] = make_IDnano_matrices_v2(Ncell, Mdisc, dt_min, gmyo, Cax, CL, CR, ggap, Gb_mat, Gc_array, length(zvec));
        [Pmax, Qmax, Cmax] = make_IDnano_matrices_v2(Ncell, Mdisc, dt_max, gmyo, Cax, CL, CR, ggap, Gb_mat, Gc_array, length(zvec));
        Pinvmin = inv(Pmin);
        Pinvmax = inv(Pmax);
end
phi = zeros(Nphi, 1);


% indexing for phi terms (in mV):
% phi_m,i:      (i-1)*(M+3)+1,          i = 1:N
% phi_R,i:      (i-1)*(M+3)+2,          i = 1:N-1
% phi_L,i:      (i-1)*(M+3),            i = 2:N
% phi_c,i,j:    (i-1)*(M+3)+j+2,        i = 1:N-1, j = 1:M
%
% indexing for I terms (in uA) and V terms (in mV): (for membrane patches)
% V/I_m,i:        (i-1)*(2*M+1)+1,        i = 1:N
% V/I_R,i,j:      (i-1)*(2*M+1)+j+1,      i = 1:N-1, j = 1:M
% V/I_L,i,j:      (i-2)*(2*M+1)+M+j+1,    i = 2:N,   j = 1:M

% scaling factor for localization of channels at ID
[~,Ncurrents] = size(loc_mat);  % number of currents in ionic model
loc_ax = 1-sum(loc_matPre,1)-sum(loc_matPost,1);  % localization on axial/lateral membrane

f_I = zeros(Ncurrents, Npatches);
for q = 1:Ncurrents
    if ~iscell(scaleI)
        for i = 1:Ncell
            f_I(q,(i-1)*(2*Mdisc+1)+1) = scaleI(q)*loc_ax(q);
            for j = 1:Mdisc
                if i>=2 % pre-junctional membrane
                    f_I(q,(i-2)*(2*Mdisc+1)+Mdisc+j+1) = scaleI(q)*loc_matPre(j,q);
                end
                
                if i<=Ncell-1  % post-junctional membrane
                    f_I(q,(i-1)*(2*Mdisc+1)+j+1) =  scaleI(q)*loc_matPost(j,q);
                end
            end
        end
    else
        if length(scaleI{q})==1
            for i = 1:Ncell
                f_I(q,(i-1)*(2*Mdisc+1)+1) = scaleI{q}*loc_ax(q);
                for j = 1:Mdisc
                    if i>=2
                        f_I(q,(i-2)*(2*Mdisc+1)+Mdisc+j+1) = scaleI{q}*loc_matPre(j,q);
                    end
                    
                    if i<=Ncell-1
                        f_I(q,(i-1)*(2*Mdisc+1)+j+1) =  scaleI{q}*loc_matPost(j,q);
                    end
                end
            end
        else
            for i = 1:Ncell
                f_I(q,(i-1)*(2*Mdisc+1)+1) = scaleI{q}(i)*loc_ax(q);
                for j = 1:Mdisc
                    if i>=2
                        f_I(q,(i-2)*(2*Mdisc+1)+Mdisc+j+1) = scaleI{q}(i)*loc_matPre(j,q);
                    end
                    
                    if i<=Ncell-1
                        f_I(q,(i-1)*(2*Mdisc+1)+j+1) =  scaleI{q}(i)*loc_matPost(j,q);
                    end
                end
            end
        end
        
    end
end
p.f_I = f_I;


G = zeros(Npatches*(p.Nstate-1),1);  % gating variables
S = zeros(3*Ncleft_comp, 1);               % cleft ion concentrations

% initial conditions
phi(((1:Ncell)-1)*(Mdisc+3)+1) = p.phii_0;  % potential
phi(((1:Ncell-1)-1)*(Mdisc+3)+2) = p.phii_0;
phi(((2:Ncell)-1)*(Mdisc+3)) = p.phii_0;

cleft_clamp = [clamp_flag(1)*ones(Ncleft_comp,1); clamp_flag(2)*ones(Ncleft_comp,1); ...
    clamp_flag(3)*ones(Ncleft_comp,1)];
Sb = [p.cleft_conc(1) p.cleft_conc(2) p.cleft_conc(3)];  % mM, bulk extracellular concentrations


G(:) = p.g0_vec; % gating variables
S(1:Ncleft_comp) = Sb(1);
S(Ncleft_comp+1:2*Ncleft_comp) = Sb(2);
S(2*Ncleft_comp+1:3*Ncleft_comp) = Sb(3);
S_bulk = S;

% indices for storing gating variables
if isempty(storeparams.indG_store) % if empty, store all of the variables
    storeparams.indG_store = 1:p.Nstate-1;
end
indGmat = [];
if ~isnan(storeparams.indG_store) % if NaN, store none of the variables
    for q = 1:length(storeparams.indG_store)
       ind_tmp = storeparams.indG_store(q);
       indGmat = [indGmat; (Npatches*(ind_tmp-1)+1:ind_tmp*Npatches)'];
    end
end

phi_mat = nan(length(phi),Ns);
Gmat = spalloc(length(G), Ns, Ns*length(indGmat));
Smat = nan(length(S), Ns);
Iall_mat = nan(Ncurrents*Npatches, Ns);
switch dtparams.method
    case 'adapt'
        ts = zeros(1, Ns);
end



t_store = 0:storeparams.int_store_all:T;
Nstore = length(t_store); count_store = 1;
phi_store = nan(length(phi), Nstore); phi_store(:,1) = phi;
G_store = nan(length(G), Nstore);  G_store(:,1) = G; 
S_store = nan(length(S), Nstore); S_store(:,1) = S;

ionic_fun = str2func(['@(t,x,p,S) ',ionic_fun_name,'(t,x,p,S)']);
Vm = constructVm(phi, Npatches, Mdisc, Ncell);

ind_Vm = ((1:Ncell)-1)*(2*Mdisc+1)+1;  % indices of axial membrane patch
ind_phim = ((1:Ncell)-1)*(Mdisc+3)+1; % indices of axial potentials
ind_phiR = ((1:Ncell-1)-1)*(Mdisc+3)+2; % indices of right potentail
ind_phiL =  ((2:Ncell)-1)*(Mdisc+3);  % indices of left potentials
ind_phic = setdiff(1:Nphi, [ind_phim ind_phiR ind_phiL]);  % indices of cleft potentials
ind_VR = []; ind_VL = [];
for i = 1:Ncell-1, for j = 1:Mdisc, ind_VR = [ind_VR, (i-1)*(2*Mdisc+1)+j+1]; end, end
for i = 2:Ncell, for j = 1:Mdisc, ind_VL = [ind_VL, (i-2)*(2*Mdisc+1)+Mdisc+j+1]; end, end


F = 96.5;                   % Faraday constant, coulombs/mmol


% intracellular diffusion
vol_ratio_ID_myo = 0.02/0.68;  % ratio of intracellular ID volume to myoplasm volume
if 0
    % from ORd11 using subspace (SS) compartment values for ID space,
    % assume 10x difference between intra ID and myoplasmic-ID diffusion
    tauMyoID_Nai = 20; tauMyoID_Cai = 20; tauMyoID_Ki = 20; % ms
    tauIntraID_Nai = 2;   tauIntraID_Cai = 2; tauIntraID_Ki = 2;
else
    ftauMyoID = 1; ftauIntraID = 5;
    L_IDID = sqrt(Ad/Mdisc); % characteristic length of intra-ID intracellular diffusion
    % diffusion parameters
    DNa = 600e-3;       % um^2 ms^-1 (from Despa and Bers 2003)
    DCa = 250e-3;       % um^2 ms^-1 (from Smith 2004)
    DK = 1300e-3;       % um^2 ms^-1 (from Hodgkin and Keynes 1953)
    tauMyoID_Nai = ftauMyoID*(L/2)^2/(2*DNa); % ms, time constant from intracellular diffusion between myoplasm and ID
    tauMyoID_Cai = ftauMyoID*(L/2)^2/(2*DCa); % ms, time constant from intracellular diffusion between myoplasm and ID
    tauMyoID_Ki = ftauMyoID*(L/2)^2/(2*DK); % ms, time constant from intracellular diffusion between myoplasm and ID
    tauIntraID_Nai = ftauIntraID*L_IDID^2/(4*DNa); % ms, time constant from intracellular diffusion between intracellular ID spaces
    tauIntraID_Cai = ftauIntraID*L_IDID^2/(4*DCa); % ms, time constant from intracellular diffusion between intracellular ID spaces
    tauIntraID_Ki = ftauIntraID*L_IDID^2/(4*DK); % ms, time constant from intracellular diffusion between intracellular ID spaces
end
% generate indices/matrices for passive intracellular diffusion coupling
[p.ind_tau_ip, p.ind_tau_im, diffusion_tau_type] = get_intracellular_diffusion_indices(Npatches, Mdisc, Gc_array, ind_Vm);
tau_Nai_mat = nan(size(p.ind_tau_ip)); tau_Cai_mat = nan(size(p.ind_tau_ip)); tau_Ki_mat = nan(size(p.ind_tau_ip));
tau_Nai_mat(diffusion_tau_type==1) = tauMyoID_Nai; tau_Nai_mat(diffusion_tau_type==2) = tauIntraID_Nai; tau_Nai_mat(diffusion_tau_type==3) = tauMyoID_Nai/vol_ratio_ID_myo;
tau_Cai_mat(diffusion_tau_type==1) = tauMyoID_Cai; tau_Cai_mat(diffusion_tau_type==2) = tauIntraID_Cai; tau_Cai_mat(diffusion_tau_type==3) = tauMyoID_Cai/vol_ratio_ID_myo;
tau_Ki_mat(diffusion_tau_type==1) = tauMyoID_Ki; tau_Ki_mat(diffusion_tau_type==2) = tauIntraID_Ki; tau_Ki_mat(diffusion_tau_type==3) = tauMyoID_Ki/vol_ratio_ID_myo;
p.tau_Nai_mat = tau_Nai_mat; p.tau_Cai_mat = tau_Cai_mat; p.tau_Ki_mat = tau_Ki_mat;


if trange(1)==0 % store initial state variables at t=0
    phi_mat(:,1) = phi; Smat(:,1) = S;
    if ~isempty(indGmat)
        Gmat(indGmat,1) = G(indGmat); 
    end
    
    Sall = nan(3*Npatches,1);
    
    Sall(ind_Vm) = Sb(1);
    Sall(ind_Vm+Npatches) = Sb(2);
    Sall(ind_Vm+2*Npatches) = Sb(3);
    for q = 1:3
        for ii = 1:Ncell-1
            Sall((q-1)*Npatches+2+(ii-1)*(2*Mdisc+1):(q-1)*Npatches+2+(ii-1)*(2*Mdisc+1)+Mdisc-1) = S((q-1)*Ncleft_comp+1+(ii-1)*Mdisc:(q-1)*Ncleft_comp+ii*Mdisc);
            Sall((q-1)*Npatches+2+(ii-1)*(2*Mdisc+1)+Mdisc:(q-1)*Npatches+ii*(2*Mdisc+1)) = S((q-1)*Ncleft_comp+1+(ii-1)*Mdisc:(q-1)*Ncleft_comp+ii*Mdisc);
        end
    end
    [~, ~, ~, ~, Iall] = ionic_fun(0,[Vm;G],p, Sall);
    Iall_mat(:,1) = Iall;
    count = 1;
else
    count = 0;
end

tic;

Vm_old = Vm;
beat_num = ones(Npatches,1);
tup = nan(Npatches,1);
trepol = nan(Npatches,1);

stimtimes = (0:floor(T/bcl))*bcl;  % ms, time of stimuli

ti = 0; H_tmp = 0; izvec = find(~clamp_flag)';
while ti < T
    
    switch dtparams.method
        case 'fixed'
            if mod(ti,bcl)<twin
                dt = dt1; dt_samp = dt1_samp; Q = Q1; C = C1; Pinv = Pinv1; dtSx = min(dt1, dtS1);
            else
                dt = dt2; dt_samp = dt2_samp; Q = Q2; C = C2; Pinv = Pinv2; dtSx = dtS2;
            end
        case {'adapt','adapt_samp'}
            if mod(ti,bcl)<twin
                dt = dt_min; Q = Qmin; C = Cmin; Pinv = Pinvmin; dtSx = min(dt_min, dtS_min); 
                if strcmp(dtparams.method,'adapt_samp')
                    dt_samp = dt_min_samp;
                end
            else
                dt = min(dt_max, dV_max/max(abs(dVdt))); dtSx = min(dtS_max, dS_max/dSdt_max); 
                if strcmp(dtparams.method,'adapt_samp')
                    dt_samp = dt_max_samp;
                end
                % determine if time step "jumps over" the next stimulus time (or end of simulation)
                ibefore = sign(stimtimes-ti);
                iafter = sign(stimtimes-(ti+dt));
                if any(ibefore.*iafter<0)
                    dt = stimtimes(ibefore.*iafter<0)-ti;
                end
                if abs(dt-dt_max)<1e-5  % if dt = dt_max
                    Q = Qmax; C = Cmax; Pinv = Pinvmax;
                else
                    [P, Q, C] = make_IDnano_matrices_v2(Ncell, Mdisc, dt, gmyo, Cax, CL, CR, ggap, Gb_mat, Gc_array, length(zvec));
                    Pinv = inv(P);
                end
            end
            
    end
    p.dt = dt; 
    Nts = ceil(dt/dtSx); dtSi = dt/Nts;  % make sure dtSi is a multiple of the current dt
   
    switch dtparams.method
        case 'fixed'
            % update  state variables
            [G_new, Iion_new, Ivec_new, Iall_new] = update_Ivec_G(ti, G, S, Vm, Npatches, Mdisc, Ncell, ...
                Ncleft_comp, Sb, p, ind_Vm,  ionic_fun);
            
            [phi_new, S_new, ~] = update_phi_S(Ivec_new, Iion_new, phi, S,  Npatches, ...
                Ncleft_comp, Sb, Pinv, Q, C, cleft_clamp, zvec, F, Vol_cleft_all,...
                Nts, dtSi,  ind_phic, ind_VL, ind_VR, Gc_array, Gb_mat, izvec);
            
        case {'adapt','adapt_samp'}
            [G_new, Iion_new, Ivec_new, Iall_new] = update_Ivec_G(ti, G, S, Vm, Npatches, Mdisc, Ncell, ...
                Ncleft_comp, Sb, p, ind_Vm,  ionic_fun);
            
            if any(isnan(G_new)) || any(abs(imag(G_new))> 1e-19) || any(isnan(Ivec_new)) || any(abs(imag(Ivec_new))> 1e-19)
                disp('Error'); break
            end

            % approximate dV
%             phi_tmp = Pinv*(Q*phi + C*Iion_new - H_tmp);
            phi_tmp = phi;  % check to make sure this is OK, defines dVdt_estimate = 0

            % determine value for dVdt
            dV = phi_tmp-phi; dVdt = dV/dt;
            
            if mod(ti,bcl)>=twin && max(abs(dV))>dV_max
                dt = max(dt_min, min(dt_max, dV_max/max(abs(dVdt))));
                ibefore = sign(stimtimes-ti);
                iafter = sign(stimtimes-(ti+dt));
                if any(ibefore.*iafter<0)
                    dt = stimtimes(ibefore.*iafter<0)-ti;
                end
                p.dt = dt;  Nts = ceil(dt/dtSx); dtSi = dt/Nts;
                if abs(dt-dt_max)<1e-5  % if dt = dt_max
                    Q = Qmax; C = Cmax; Pinv = Pinvmax;
                else
                    [P, Q, C] = make_IDnano_matrices_v2(Ncell, Mdisc, dt, gmyo, Cax, CL, CR, ggap, Gb_mat, Gc_array, length(zvec));
                    Pinv = inv(P);
                end
                % update  state variables with new dt
                [G_new, Iion_new, Ivec_new, Iall_new] = update_Ivec_G(ti, G, S, Vm, Npatches, Mdisc, Ncell, ...
                    Ncleft_comp, Sb, p, ind_Vm, ionic_fun);
                
                if any(isnan(G_new)) || any(abs(imag(G_new))> 1e-19) || any(isnan(Ivec_new)) || any(abs(imag(Ivec_new))> 1e-19)
                    disp('Error'); break;
                end

                
                [phi_new, S_new, H_tmp, dSdt_max] = update_phi_S(Ivec_new, Iion_new, phi, S,  Npatches, ...
                    Ncleft_comp, Sb, Pinv, Q, C, cleft_clamp, zvec, F, Vol_cleft_all,...
                    Nts, dtSi, ind_phic, ind_VL, ind_VR, Gc_array, Gb_mat, izvec);
                            
            else
                [phi_new, S_new, H_tmp, dSdt_max] = update_phi_S(Ivec_new, Iion_new, phi, S,  Npatches, ...
                    Ncleft_comp, Sb, Pinv, Q, C, cleft_clamp, zvec, F, Vol_cleft_all,...
                    Nts, dtSi,  ind_phic, ind_VL, ind_VR, Gc_array, Gb_mat, izvec);
            end
            % determine value for dVdt
            dV = phi_new-phi; dVdt = dV/dt; % for next iteration
    end

    
    % construct vector of potential differences, Vm_new
    Vm = constructVm(phi_new, Npatches, Mdisc, Ncell);
    
    if any(isnan(phi_new)) || any(abs(imag(phi_new))> 1e-19)
        disp('Error'); break;
    end
    
    % calculate activation/repolarization times
    [tup, trepol, beat_num] = update_tup_repol(ti, dt, Vm, Vm_old, Vthresh, ...
        tup, trepol, beat_num);
    
    if single_beat_flag && ~any(isnan(tup(:,1)))
        toc; 
        tup_Vm = tup(ind_Vm,:);
        i1 = round(.25*Ncell); i2 = round(.75*Ncell);
        cv_est = 100*(i2-i1)*(p.L/1000)./(tup_Vm(i2,:)-tup_Vm(i1,:));  % mm/ms = m/s, 100*m/s = cm/s
        
        store_out.phi_store = phi_store;
        store_out.G_store = G_store;
        store_out.S_store = S_store;
        store_out.t_store = t_store;
        return;
    end
    Vm_old = Vm;
    
    ti = round(ti + dt, 5);
    
    % save values
    switch dtparams.method
        case 'fixed'
            if ti>=trange(1) && ti<=trange(2) && ~mod(ti, dt_samp)
                count = count + 1;
                phi_mat(:,count) = phi_new;
                if ~isempty(indGmat)
                    Gmat(indGmat,count) = G_new(indGmat);
                end
                Smat(:,count) = S_new;
                Iall_mat(:,count) = Iall_new;
            end
        case 'adapt'
            if ti>=trange(1) && ti<=trange(2)
                count = count + 1;
                ts(count) = ti;
                phi_mat(:,count) = phi_new;
                if ~isempty(indGmat)
                    Gmat(indGmat,count) = G_new(indGmat);
                end
                Smat(:,count) = S_new;
                Iall_mat(:,count) = Iall_new;
            end
        case 'adapt_samp'
            if ti>=trange(1) && ti<=trange(2)
                if ~mod(ti, dt_samp)
                    count = count + 1;
                    phi_mat(:,count) = phi_new;
                    if ~isempty(indGmat)
                        Gmat(indGmat,count) = G_new(indGmat);
                    end
                    Smat(:,count) = S_new;
                    Iall_mat(:,count) = Iall_new;
                else
                    ibefore = sign(ts-(round(ti-dt,5)));
                    iafter = sign(ts-ti);
                    if any(ibefore.*iafter<0)
                        count = count + 1;
                        % linear interpolate between old and new data points
                        tx = ts(ibefore.*iafter<0);
                        
                        mphi = (phi_new-phi)/dt;
                        phi_mat(:,count) = phi_new + mphi*(tx-ti);
                        
                        mG = (G_new-G)/dt;
                        if ~isempty(indGmat)
                            Gmat(indGmat,count) = G_new(indGmat) + mG(indGmat)*(tx-ti);
                        end
                        mS = (S_new-S)/dt;
                        Smat(:,count) = S_new + mS*(tx-ti);
                        
                        mI = (Iall_new-Iall)/dt;
                        Iall_mat(:,count) = Iall_new + mI*(tx-ti);
                        
                    end
                end  
            end
    end
    
    % store values
    if ~mod(ti, storeparams.int_store_all)
        phi_store(:,count_store+1) = phi_new;
        G_store(:,count_store+1) = G_new;
        S_store(:,count_store+1) = S_new;
        count_store = count_store + 1;

    end
        
    G = G_new;
    phi = phi_new;
    S = S_new;
    Iall = Iall_new;
    
    % plot during simulation
        if 0 && ~mod(count-1,50)
            figure(100);
    
            phi_i = phi_mat(ind_phim,:);
            phi_cleft = phi_mat(ind_phic,:);
            subplot(2,2,1); plot(ts, phi_i);
            subplot(2,2,2); plot(ts, phi_cleft);
            subplot(2,2,3); plot(ts, Smat(1:Ncleft_comp,:));
            subplot(2,2,4); plot(ts, Smat(Ncleft_comp + (1:Ncleft_comp),:));
           % subplot(2,3,6); plot(ts, Smat(2*Ncleft_comp + (1:Ncleft_comp),:));
            %         imagesc(1:Ncell, ts, phi_i'-phi_e');
            %         imagesc(1:Ncell, ts, phi_i');
            %         title(['t = ',num2str(t(i))]);
            %         set(gca,'ydir','normal');
            drawnow;
        end
end

switch dtparams.method
    case 'adapt'
        ts = ts(1:count);
        phi_mat = phi_mat(:,1:count);
        Gmat = Gmat(:,1:count);
        Smat = Smat(:,1:count);
end

toc;

p.Ncurrents = Ncurrents;
p.Npatches = Npatches;
p.Ns = Ns;

store_out.phi_store = phi_store;
store_out.G_store = G_store;
store_out.S_store = S_store;
store_out.t_store = t_store;

if nargin == 0
    phi_i = phi_mat(ind_phim,:);
  %  phi_cleft = phi_mat(ind_phic,:);
    
    figure; imagesc(1:Ncell, ts, phi_i');  set(gca,'ydir','normal');
    
end

%  conduction velocity calculation (deltaX/deltaT)
% use middle 50% of cable for CV estimate
tup_Vm = tup(ind_Vm,:);
i1 = round(.25*Ncell); i2 = round(.75*Ncell);    
cv_est = 100*(i2-i1)*(p.L/1000)./(tup_Vm(i2,:)-tup_Vm(i1,:));  % mm/ms = m/s, 100*m/s = cm/s


end


function [tup, trepol, beat_num] = update_tup_repol(ti, dt, Vm, Vm_old, Vthresh, ...
    tup, trepol, beat_num)
iact = find((Vm>Vthresh).*(Vm_old<Vthresh));
for j = 1:length(iact)
    y1 = Vm_old(iact(j)); y2 = Vm(iact(j));
    m = (y2-y1)/dt;
    tup(iact(j),beat_num(iact(j))) = ti - (y1-Vthresh)./m;
    %         save(times_name,'tup','trepol');
end
irepol = find((Vm<Vthresh).*(Vm_old>Vthresh));
for j = 1:length(irepol)
    if ti-tup(irepol(j),beat_num(irepol(j)))>.1
        y1 = Vm_old(irepol(j)); y2 = Vm(irepol(j));
        m = (y2-y1)/dt;
        trepol(irepol(j),beat_num(irepol(j))) = ti - (y1-Vthresh)./m;
        beat_num(irepol(j)) = beat_num(irepol(j)) + 1;
    end
    %         save(times_name,'tup','trepol');
end

end


function Vm = constructVm(phi, Npatches, Mdisc, Ncell)

% construct vector of potential differences
Vm = nan(Npatches,1);
% indexing for phi terms (in mV):
% phi_m,i:      (i-1)*(M+3)+1,          i = 1:N
% phi_R,i:      (i-1)*(M+3)+2,          i = 1:N-1
% phi_L,i:      (i-1)*(M+3),            i = 2:N
% phi_c,i,j:    (i-1)*(M+3)+j+2,        i = 1:N-1, j = 1:M

% indexing for I terms (in uA) and V terms (in mV): (for membrane patches)
% V/I_m,i:        (i-1)*(2*M+1)+1,        i = 1:N
% V/I_R,i,j:      (i-1)*(2*M+1)+j+1,      i = 1:N-1, j = 1:M
% V/I_L,i,j:      (i-2)*(2*M+1)+M+j+1,    i = 2:N,   j = 1:M

ind_Vm = ((1:Ncell)-1)*(2*Mdisc+1)+1;
ind_phim = ((1:Ncell)-1)*(Mdisc+3)+1;
for i = 1:Ncell
    
    if i<=Ncell-1
        ind_phiR = (i-1)*(Mdisc+3)+2;
        for j = 1:Mdisc
            ind_phic = (i-1)*(Mdisc+3)+j+2;
            ind_VR = (i-1)*(2*Mdisc+1)+j+1;
            
            Vm(ind_VR) = phi(ind_phiR) - phi(ind_phic);  % VR
        end
    end
    
    if i>=2
        ind_phiL = (i-1)*(Mdisc+3);
        for j = 1:Mdisc
            ind_phic = (i-1-1)*(Mdisc+3)+j+2;
            ind_VL = (i-2)*(2*Mdisc+1)+Mdisc+j+1;
            
            Vm(ind_VL) = phi(ind_phiL) - phi(ind_phic);  % VL
        end
    end
end
Vm(ind_Vm) = phi(ind_phim);  % Vm_i

end

function [ind_tau_ip, ind_tau_im, diffusion_tau_type] = get_intracellular_diffusion_indices(Npatches, Mdisc, Gc_array, ind_Vm)

ind_tau_ip = nan(Npatches, 2*Mdisc);
ind_tau_im = nan(Npatches, 2*Mdisc);
diffusion_tau_type = nan(Npatches, 2*Mdisc);  % 1 = (Myo->ID), 2 = (ID<->ID), 3 = (Myo<-ID)
iGci = sign(Gc_array(:,:,1));
for i = 1:Npatches
    ind_tau_im(i,:) = i;
    ind_tau_ip(i,:) = i;
end
for i = 1:length(ind_Vm)
    ind_tau_ip(ind_Vm(i),1:Mdisc) = max(1,ind_Vm(i)-Mdisc:ind_Vm(i)-1);
    ind_tau_ip(ind_Vm(i),1+Mdisc:2*Mdisc) = min(Npatches,ind_Vm(i)+1:ind_Vm(i)+Mdisc);
    diffusion_tau_type(ind_Vm(i),:) = 3;
    
    if ind_Vm(i) == 1
        ind = ind_Vm(i)+1:ind_Vm(i)+Mdisc;
        for j = 1:Mdisc
            ind_j = ind(j);
            icoup = find(iGci(j,:));
            ind_tau_ip(ind_j,1:length(icoup)) = ind(icoup);
            diffusion_tau_type(ind_j,:) = 2;
            ind_tau_ip(ind_j,length(icoup)+1) = ind_Vm(i);
            diffusion_tau_type(ind_j, length(icoup)+1) = 1;
        end
    elseif ind_Vm(i) == Npatches
        ind = ind_Vm(i)-Mdisc:ind_Vm(i)-1;
        for j = 1:Mdisc
            ind_j = ind(j);
            icoup = find(iGci(j,:));
            ind_tau_ip(ind_j,1:length(icoup)) = ind(icoup);
            diffusion_tau_type(ind_j,:) = 2;
            ind_tau_ip(ind_j,length(icoup)+1) = ind_Vm(i);
            diffusion_tau_type(ind_j, length(icoup)+1) = 1;
        end
    else
        % left side
        ind = ind_Vm(i)-Mdisc:ind_Vm(i)-1;
        for j = 1:Mdisc
            ind_j = ind(j);
            icoup = find(iGci(j,:));
            ind_tau_ip(ind_j,1:length(icoup)) = ind(icoup);
            diffusion_tau_type(ind_j,:) = 2;
            ind_tau_ip(ind_j,length(icoup)+1) = ind_Vm(i);
            diffusion_tau_type(ind_j, length(icoup)+1) = 1;
        end
        % right side
        ind = ind_Vm(i)+1:ind_Vm(i)+Mdisc;
        for j = 1:Mdisc
            ind_j = ind(j);
            icoup = find(iGci(j,:));
            ind_tau_ip(ind_j,1:length(icoup)) = ind(icoup);
            diffusion_tau_type(ind_j,:) = 2;
            ind_tau_ip(ind_j,length(icoup)+1) = ind_Vm(i);
            diffusion_tau_type(ind_j, length(icoup)+1) = 1;
        end
    end
end
end

function ts = get_time_variable(trange, dt1, dt2, twin, bcl)

ti = trange(1);  count = 1;
ts = nan(1,1e6);
ts(1) = trange(1);
while ti < trange(2)
    if mod(ti,bcl)<twin
        dt = dt1;
    else
        dt = dt2;
    end
    
    ti = round(ti + dt,4);
    ts(count+1) = ti;
    count = count + 1;
end
ts = ts(1:count);
end

function [G, Iion, Ivec, Iall] = update_Ivec_G(ti, G, S, Vm, Npatches, Mdisc, Ncell, ...
    Ncleft_comp, Sb, p, ind_Vm, ionic_fun)

% construct extracellular ion concentration vector (assume lateral
% membrane extracellular ion concentrations are constant)
Sall = nan(3*Npatches,1);

Sall(ind_Vm) = Sb(1);
Sall(ind_Vm+Npatches) = Sb(2);
Sall(ind_Vm+2*Npatches) = Sb(3);
for q = 1:3
    for ii = 1:Ncell-1
        Sall((q-1)*Npatches+2+(ii-1)*(2*Mdisc+1):(q-1)*Npatches+2+(ii-1)*(2*Mdisc+1)+Mdisc-1) = S((q-1)*Ncleft_comp+1+(ii-1)*Mdisc:(q-1)*Ncleft_comp+ii*Mdisc);
        Sall((q-1)*Npatches+2+(ii-1)*(2*Mdisc+1)+Mdisc:(q-1)*Npatches+ii*(2*Mdisc+1)) = S((q-1)*Ncleft_comp+1+(ii-1)*Mdisc:(q-1)*Ncleft_comp+ii*Mdisc);
    end
end

% update G_new, Iion_new -> calc using G_old, phi_old (in Vm terms), S_old
[G, Iion, Ivec, ~, Iall] = ionic_fun(ti,[Vm;G],p, Sall);  % update ionic currents, gating variables

end

function [phi, S, H, dSdt_max] = update_phi_S(Ivec, Iion, phi, S,  Npatches, ...
    Ncleft_comp, Sb, Pinv, Q, C, cleft_clamp, zvec, F, Vol_cleft_all,...
    Nts, dtSi,  ind_phic, ind_VL, ind_VR, Gc_array, Gb_mat, izvec)

dS = zeros(size(S));
% update S_new -> calc using phi_old, S_old, Iion_new (actually Iion_old)
for qq = 1:Nts
    [H, Icleft] = calc_cleft_currents(S, phi(ind_phic), Gc_array, Gb_mat, Sb, zvec, izvec); % update cleft currents
    
    for q = izvec % loop over all ionic species
        ind_S = (q-1)*Ncleft_comp + (1:Ncleft_comp);
        ind_patch = (q-1)*Npatches + (1:Npatches);
        ind_L = ind_patch(ind_VL);
        ind_R = ind_patch(ind_VR);
        
        Iall_q = Ivec(ind_L) + Ivec(ind_R) - Icleft(ind_S);
        dS(ind_S) = ~cleft_clamp(ind_S)*dtSi.*(Iall_q*1e6./(zvec(q)*F*Vol_cleft_all)); % convert uA/um^3 to mM/msec
        S(ind_S) = S(ind_S) + dS(ind_S);
            
    end
    
    if any(isnan(S)) || any(abs(imag(S))> 1e-19) || any(S<0)
        disp('Error'); 
    end
    
end

dSdt_max = max(abs(dS))/dtSi;

% update phi_new -> calc using phi_old, S_new (in H term), Iion_new
%         phi = P\(Q*phi + C*Iion - H);
phi = Pinv*(Q*phi + C*Iion - H);

end