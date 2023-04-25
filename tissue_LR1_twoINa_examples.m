
Ncell = 50;     % number of cells
% cell geometry parameters
L = 100;        % cell length, um
r = 11;         % cell radius, um

frac2 = .7;  % fNa
mix_frac = [1-frac2 frac2];
total_fNa = 1; % total INa conductance scaling factor

w = 10e-3; % cleft width, um
Dscale = 1;  % 1 corresponds with 100 nS, for r = 11 um, L = 100 um
D_val = Dscale*.131; % .131

bcl = 500; % basic cycle length, ms
nbeats = 2; % number of beats

D = D_val*ones(Ncell-1,1); % diffusion coefficient, cm^2/s
ggap = convert_D_to_ggap(D_val, r, L); % gap junction conducance, mS (based on geometry)

Ncell_vec0 = [.75 1 1.25];
Ncell_vec = Ncell_vec0;
for i = 2:Ncell
    Ncell_vec = [Ncell_vec, Ncell_vec(end-2:end)+ 1];
end
Ncell_vec = Ncell_vec(2:end-1);
[~,ind1] = min(abs(Ncell_vec-24.25));
[~,ind2] = min(abs(Ncell_vec-24.75));
[~,ind3] = min(abs(Ncell_vec-25));

single_beat_flag = 0;

s1 = 'kb';
s2 = 'gr';
s = 'km';
titl = {'Baseline/Polarized','Mix/Het-Polarized'};
for q = 1:2
    switch q
        case 1 % baseline, polarized
            fNa1 = 1*total_fNa; % INa scaling factor (baseline)
            fNa2 = 0; % INa scaling factor (shifted)
            locINa1 = mix_frac(2);  % overall localization at ID
            locINa2 = 0;  % overall localization at ID

        case 2 % mix, hom-polarized
            fNa1 = total_fNa*mix_frac(1); % INa scaling factor (baseline)
            fNa2 = total_fNa*mix_frac(2); % INa scaling factor (shifted)
            locINa1 = 0;  % overall localization at ID
            locINa2 = 1;  % overall localization at ID

        case 3 % all baseline Na channels, all lateral
            fNa1 = 1*total_fNa; % INa scaling factor (baseline)
            fNa2 = 0; % INa scaling factor (shifted)
            locINa1 = 0;  % overall localization at ID
            locINa2 = 0;  % overall localization at ID

        case 4 % all shifted Na channels, all lateral
            fNa1 = 0; % INa scaling factor (baseline)
            fNa2 = 1*total_fNa; % INa scaling factor (shifted)
            locINa1 = 0;  % overall localization at ID
            locINa2 = 0;  % overall localization at ID

        case 5 % mix of channels, all lateral
            fNa1 = total_fNa*mix_frac(1); % INa scaling factor (baseline)
            fNa2 = total_fNa*mix_frac(2); % INa scaling factor (shifted)
            locINa1 = 0;  % overall localization at ID
            locINa2 = 0;  % overall localization at ID

        case 6 % all shifted channels, polarized
            fNa1 = 0; % INa scaling factor (baseline)
            fNa2 = 1*total_fNa; % INa scaling factor (shifted)
            locINa1 = 0;  % overall localization at ID
            locINa2 = mix_frac(2);  % overall localization at ID

        case 7 % mix, het-polarized (heterogeneous)
            fNa1 = total_fNa*mix_frac(1); % INa scaling factor (baseline)
            fNa2 = total_fNa*mix_frac(2); % INa scaling factor (shifted)
            locINa1 = mix_frac(2);  % overall localization at ID (note: this is not a typo, should be mix_frac(2))
            locINa2 = mix_frac(2);  % overall localization at ID


    end

    [cv, tup, ts, phi_i, phi_cleft, phi_L, phi_R, INa1, INa2]...
        = fun_study_run_LR1_twoINa_epc_singecleft(Ncell, D, fNa1, fNa2, ...
        w, locINa1, locINa2, L, r, total_fNa, single_beat_flag, bcl, nbeats);

    cv

    if Dscale < 1
        tup = tup(:,1);
    end

    figure(1);
    subplot(3,3,1);
    plot(Ncell_vec, tup(:,end),'-','color',s(q),'linewidth',2,'markersize',18); hold on;

    subplot(3,3,4);
    plot(Ncell_vec, tup(:,end)-tup(ind1,end),'.-','color',s(q),'linewidth',2,'markersize',18); hold on;
    xlim([24 26]);

    subplot(3,3,7);
    plot(Ncell_vec, tup(:,end)-tup(ind2,end),'.-','color',s(q),'linewidth',2,'markersize',18); hold on;
    xlim([24 26]);

    icell = 25;
    subplot(4,3,q+1);
    plot(ts,phi_i(icell-1:icell+1,:),s1(q),'linewidth',1.5); hold on;
    plot(ts,phi_R(icell-1:icell+1,:),[s2(q),'--'],'linewidth',1.5)
    plot(ts,phi_L(icell-2:icell,:),[s2(q),'--'],'linewidth',1.5);
    xlim(tup(ind3,end)+[-.5 1]); ylim([-100 30])
    plot(tup(ind3,end)+[.4 .9],-95*[1 1],'k','linewidth',8)
    text(tup(ind3,end)+.65,-108,'0.5 ms','HorizontalAlignment','center','fontsize',24);
    set(gca,'fontsize',28,'tickdir','out','linewidth',2,'box','off');
    if q == 1
        ylabel('\phi (mV)');
    else
        set(gca,'ytick',[],'ycolor','w')
    end
    title(titl{q});
    set(gca,'xtick',[],'xcolor','w')

    subplot(4,3,q+4);
    plot(ts, phi_cleft(icell-1:icell,:),s(q),'linewidth',2); hold on;
    set(gca,'fontsize',28,'tickdir','out','linewidth',2,'box','off');
    xlim(tup(ind3,end)+[-.5 1]); ylim([-30 2]);
    set(gca,'xtick',[],'xcolor','w')
    if q == 1
        ylabel('\phi_{cleft} (mV)');
    else
        set(gca,'ytick',[],'ycolor','w')
    end

    subplot(4,3,q+7);
    switch q
        case 1 %  baseline, polarized
            plot(ts, 1e3*INa1([ind1-1 ind1+2 ind1+5],:),'k','linewidth',1.5); hold on;
            plot(ts, 1e3*INa1([ind1 ind1+3 ],:),'g','linewidth',1.5); hold on;
            plot(ts, 1e3*INa1([ind1+1  ind1+4],:),'g--','linewidth',1.5); hold on;

        case 2 % mix, hom-polarized
            plot(ts, 1e3*INa1([ind1-1 ind1+2 ind1+5],:),'b','linewidth',1.5); hold on;
            plot(ts, 1e3*INa2([ind1 ind1+3 ],:),'r','linewidth',1.5); hold on;
            plot(ts, 1e3*INa2([ind1+1  ind1+4],:),'r--','linewidth',1.5); hold on;

    end
    set(gca,'fontsize',28,'tickdir','out','linewidth',2,'box','off');
    xlim(tup(ind3,end)+[-.5 1]);
    ylim([-.015*1e3 0]);
    set(gca,'xtick',[],'xcolor','w')
    if q == 1
        ylabel('I_{Na} (nA)');
    else
        set(gca,'ytick',[],'ycolor','w')
    end

    subplot(4,3,q+10);
    Igap = 1e3*ggap*(phi_R(icell-1:icell,:)-phi_L(icell-1:icell,:)); % uA
    plot(ts, Igap,s(q),'linewidth',2); hold on;
    set(gca,'fontsize',28,'tickdir','out','linewidth',2,'box','off');
    xlim(tup(ind3,end)+[-.5 1]); ylim([-1 8]);
    set(gca,'xtick',[],'xcolor','w')
    if q == 1
        ylabel('I_{gap} (nA)');
    else
        set(gca,'ytick',[],'ycolor','w')
    end


end

figure(1)
subplot(3,3,1); xlim([0 50.5])
set(gca,'fontsize',28,'tickdir','out','linewidth',2,'box','off');
ylabel('Activation time (ms)'); %xlabel('Cell #');
h = legend('Baseline/Polarized','Mix/Hom-Polarized');
h.EdgeColor = 'w';

subplot(3,3,4);
plot(24.5*[1 1],[-1 1],'k--','linewidth',1); ylim([-0.02 1])
plot(25.5*[1 1],[-1 1],'k--','linewidth',1);
set(gca,'fontsize',28,'tickdir','out','linewidth',2,'box','off');
ylabel('\Delta Activation time (ms)'); %xlabel('Cell #');


subplot(3,3,7);
plot(24.5*[1 1],[-1 1],'k--','linewidth',1); ylim([-0.22 .8])
plot(25.5*[1 1],[-1 1],'k--','linewidth',1);
set(gca,'fontsize',28,'tickdir','out','linewidth',2,'box','off');
ylabel('\Delta Activation time (ms)'); xlabel('Cell #');

