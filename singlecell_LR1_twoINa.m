
% Units
%
% voltage in mV
% current in uA
% conductance in mS
% resistance in k-ohm
% capacitance in uF    (uF*mV/ms = uA)
% time in ms
% length in um
% concentration in mM
%


s = 'k-'; clear p;
model = 'LR1_twoINa';

fINa1 = .3;
fINa2 = .7;


for qqq = [2 1]
    switch qqq
        case 1
            Vact_shift = 5.5; % mV, shift in SS activation curve (+ means left shift)
            Vinact_shift  = -7.8; % mV, shift in SS inactivation curve, (- means right)
        case 2
            Vact_shift = 0;
            Vinact_shift = 0;
    end

    % cell geometry parameters
    L = 100;        % cell length, um
    r = 11;         % cell radius, um


    % time parameters
    bcl = 500;      % basic cycle length, ms
    nbeats = 10;      % number of beats
    T = bcl*nbeats;   % total time, ms

    % cell geometry
    Aax = 2*pi*r*L; % axial patch surface area, um^2
    Ad = pi*r^2;    % disc surface area, um^2
    Atot = 2*Ad + Aax;  % total surface area, um^2

    % bulk extracellular concentrations
    Ko = 5.4;                  % mM
    Nao = 140;                 % mM
    Cao = 1.8;                 % mM

    Npatches = 1; % number of membrane patches

    % model specific parameters
    switch model

        case 'LR1_twoINa'
            % initial conditions
            [p, x0] = InitialConstants_LR91_twoINa(Atot);

            % order is determined by code in fun_name
            % INa1, Isi, IK, IK1, IKp, Ib, INa2
            Ncurrents = 7;
            scaleI = ones(1,Ncurrents);

            % ionic model-specific parameters
            ionic_fun_name = 'fun_LR1_twoINa';

            p.iina1 = 1; p.iisi = 2; p.iik = 3;
            p.iik1 = 4;  p.iikp = 5; p.iib = 6;
            p.iina2 = 7;

            %             p.Vact_shift = 5.5; % mV, shift in SS activation curve (+ means left shift)
            %             p.Vinact_shift  = -7.8; % mV, shift in SS inactivation curve, (- means right)

            p.Vact_shift = Vact_shift;
            p.Vinact_shift = Vinact_shift;
            scaleI(p.iina1) = fINa1;  % example of doubling the INa conductance
            scaleI(p.iina2) = fINa2;
            p.Nstate = 8+3;  % number of state variables, including Vm, per patch

            % model specific parameter
            p.mLR1 = 1; % flag for modified LR1 model with Ca2+ speedup

            % stimulus parameters
            p.stim_dur = 1;   % ms
            p.stim_amp = 80e-8*Atot;    % uA

    end

    p.L = L; p.r = r;  % um
    % extracellular concentrations
    p.K_o = Ko;                  % mM
    p.Na_o = Nao;                 % mM
    p.Ca_o = Cao;                 % mM
    % cleft ionic concentration parameters
    p.cleft_conc = [p.Na_o; p.K_o; p.Ca_o];
    zvec = [1; 1; 2];  % charge valence

    psingle = p;
    psingle.bcl = bcl;
    psingle.Npatches = 1;
    psingle.f_I = scaleI';

    psingle.indstim = 1;
    psingle.celltype = 0;
    psingle.ind_tau_ip = 1; psingle.ind_tau_im = 1; psingle.tau_Nai_mat = 1; psingle.tau_Cai_mat = 1; psingle.tau_Ki_mat = 1;

    Cm = 1*1e-8;      % membrane capacitance, uF/um^2
    psingle.Ctot = (2*Ad + 2*pi*r*L)*Cm;   % total cell capacitance, uF


    Vm_single(1) = x0(1); % assumes the first variable is always Vm
    Gsingle = x0(2:end)';

    dt = 1e-2; % time step, ms
    psingle.dt = dt;

    options = odeset('MaxStep',1,'InitialStep',2e-2);

    opt = odeset('maxstep',1);
    tic;
    ionic_fun = str2func(['@(t,x,p,S) ',ionic_fun_name,'(t,x,p,S)']);
    [t,X] = ode15s(@(t,x) get_ionic_fun_dX(ionic_fun,t,x,psingle, p.cleft_conc),...
        0:dt:bcl*nbeats,[Vm_single;Gsingle], options);

   
    toc;

    % calculate ionic currents
    Ncurrents = length(scaleI);
    Iall_mat = nan(Ncurrents, length(t));
    for i = 1:length(t)
        [~, ~, ~, ~, Iall_mat(:,i)] = ionic_fun(t(i), X(i,:), psingle, p.cleft_conc);
    end

    [~,ind] = min(abs(t-(nbeats-1)*bcl));
    ratio = min(Iall_mat(p.iina2,ind:end))./min(Iall_mat(p.iina1,ind:end));

    s = 'mk';
    V = X(:,1);
    subplot(3,1,1); plot(t-(nbeats-1)*bcl, V,s(qqq),'linewidth',2); hold on; xlim([-10 400])
    title([num2str(100*fINa1),'-',num2str(100*fINa2),' I_{Na} distribution'])

    subplot(3,1,2); plot(t-(nbeats-1)*bcl, V,s(qqq),'linewidth',2); hold on; xlim([-.1 1])

    switch qqq
        case 1
            subplot(3,1,3);
            plot(t-(nbeats-1)*bcl, 1e3*Iall_mat(p.iina1,:),'b','linewidth',2); hold on;
            plot(t-(nbeats-1)*bcl, 1e3*Iall_mat(p.iina2,:),'r','linewidth',2); hold on;
            plot(t-(nbeats-1)*bcl, 1e3*(Iall_mat(p.iina2,:)+Iall_mat(p.iina1,:)),'m','linewidth',2); hold on;
            xlim([-.1 1])
        case 2
            subplot(3,1,3);
            plot(t-(nbeats-1)*bcl, 1e3*(Iall_mat(p.iina1,:)+Iall_mat(p.iina2,:)),'k','linewidth',2); hold on;
            
    end
end

subplot(3,1,1); set(gca,'fontsize',28,'tickdir','out','linewidth',2,'box','off')
ylabel('V (mV)'); xlabel('time (ms)'); ylim([-100 80])

subplot(3,1,2); set(gca,'fontsize',28,'tickdir','out','linewidth',2,'box','off','xcolor','w')
ylabel('V (mV)'); %ylabel('time (ms)')
h = legend('Homogeneous (Only Baseline)','Baseline/Shifted'); h.EdgeColor = 'w'; ylim([-100 80])

subplot(3,1,3); set(gca,'fontsize',28,'tickdir','out','linewidth',2,'box','off')
ylabel('I_{Na} (nA)'); xlabel('time (ms)');
ylim([-35 0])
h = legend('Homogenous I_{Na}','Baseline I_{Na}','Shifted I_{Na}','Total Baseline/Shifted I_{Na}');
h.EdgeColor = 'w';


