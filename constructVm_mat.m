function Vm = constructVm_mat(phi, Npatches, Mdisc, Ncell)
[~,nt] = size(phi);
% construct vector of potential differences
Vm = nan(Npatches,nt);
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
            
            Vm(ind_VR,:) = phi(ind_phiR,:) - phi(ind_phic,:);  % VR
        end
    end
    
    if i>=2
        ind_phiL = (i-1)*(Mdisc+3);
        for j = 1:Mdisc
            ind_phic = (i-1-1)*(Mdisc+3)+j+2;
            ind_VL = (i-2)*(2*Mdisc+1)+Mdisc+j+1;
            
            Vm(ind_VL,:) = phi(ind_phiL,:) - phi(ind_phic,:);  % VL
        end
    end
end
Vm(ind_Vm,:) = phi(ind_phim,:);  % Vm_i

end