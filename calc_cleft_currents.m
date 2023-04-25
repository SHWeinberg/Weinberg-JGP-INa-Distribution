function [H, Icleft] = calc_cleft_currents(S, phic, Gc_array, Gb_mat, Sb, zvec, izvec)

% inputs:
% S: mM, cleft ion concentrations, 3*Ncomp_cleft x 1 vector, Ncomp_cleft = (Ncell-1)*Mdisc
% phic: mV, cleft potentials, Ncomp_cleft x 1 vector
% Gc_array: mS, cleft conductances, Mdisc x Mdisc x Ncell array, Gc(j,k,i) = conductance between patch j and k for cleft i
% Gb_mat: mS, bulk conductance, Ncell x Mdisc matrix
% Sb: mM, bulk ion concentrations, 3 x 1 vector (Na+, K+, Ca2+)
% zvec: unitless, valence vector, 3 x 1
%
% outputs:
% H: uA, summation of reversal potential terms, Nphi x 1 vector
% Icleft: uA, summation of individual ionic species currents between cleft
% compartments, 3*Ncomp_cleft x 1 (to be used charge conservation equation)


[Ncell, Mdisc] = size(Gb_mat);
Ncleft_comp = length(phic);
Nphi = 3*Ncell-2 + Mdisc*(Ncell-1);  % number of voltage state variables

F = 96.5;                   % Faraday constant, coulombs/mmol
R = 8.314;                  % gas constant, J/K
T = 273+37;                 % absolute temperature, K
RTF=(R*T/F);                % mV

H_z = zeros(Nphi, 3);
Icleft = zeros(length(zvec)*Ncleft_comp, 1);


for i = 1:Ncell-1
    for j = 1:Mdisc
        r = j + (i-1)*Mdisc + 3*Ncell - 2;  % EQN. 4
        ind_j = (i-1)*Mdisc + j;

        for q = izvec  % izvec, loop over each ionic species (hardcoding speeds up by ~10% for more than 1 species)
            ind_j_q = (q-1)*Ncleft_comp + (i-1)*Mdisc + j;
            Erevj = RTF/zvec(q)*log(Sb(q)/S(ind_j_q));
            
            gc_Eterms = zeros(1,Mdisc);
            gc_Iterms = zeros(1,Mdisc);
                
            for k = find(Gc_array(j,:,i))
                ind_k_q = (q-1)*Ncleft_comp + (i-1)*Mdisc + k;
                ind_k = (i-1)*Mdisc + k;
                
                Erevkj = RTF/zvec(q)*log(S(ind_k_q)/S(ind_j_q));
                
                gc_Eterms(k) = Gc_array(j,k,i)*Erevkj;
                gc_Iterms(k) = Gc_array(j,k,i)*(phic(ind_j) - phic(ind_k) - Erevkj);
            end
            H_z(r,q) = sum(gc_Eterms) + Gb_mat(i,j)*Erevj;
            Icleft(ind_j + (q-1)*Ncleft_comp) = sum(gc_Iterms) + Gb_mat(i,j)*(phic(ind_j) - Erevj);
        end
    end
    
end
H = sum(H_z,2);