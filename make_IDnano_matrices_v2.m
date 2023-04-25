function [P, Q, C] = make_IDnano_matrices_v2(N, M, dt, gmyo, Cm, CL, CR, ggap, Gb, Gc, ns)

if nargin == 0
    M = 20;      % number of cleft patches
    N = 50;      % number of cells
    
    % Units
    %
    % voltage in mV
    % current in uA  (mS*mV = uA)
    % conductance in mS, resistance in k-ohm
    % capacitance in uF    (uF/ms = mS)
    % time in ms
    
    dt = 1;  
    
    gmyo = 2; % mS, myoplasmic conductance, scalar
    Cm = 3;   % uF, axial membrane capacitance, scalar
    
    CL = 5*ones(N,M);  % uF, left side capacitance, N x M matrix
    CR = 7*ones(N,M);  % uF, right side capacitance, N x M matrix
    
    ggap = 11*ones(N-1,1); % mS, gap junctional conductance, N-1 x 1 vector, between cell i and i+1, for i = 1,...,N-1
    
    Gb = 13*ones(N, M); % mS, bulk conductance, N x M matrix
    % for cleft i (right of cell i), Gb(i,j) = conductance between patch j and bulk
    
    Gci = 17*ones(M,M);
    for j = 1:M, Gci(j,j) = 0; end
    Gc = nan(M,M,N);
    % mS, cleft conductances, M x M x N array, Gc(j,k,i) = conductance between patch j and k for cleft i
    for i = 1:N, Gc(:,:,i) = Gci; end
    
    ns = 3; % number of ionic species, i.e., Na+, K+, Ca2+
end

Nphi = 3*N-2 + M*(N-1);
Npatches = 2*M*(N-1) + N;
Ncleft_comp = M*(N-1);

P = sparse(Nphi, Nphi);
Q = sparse(Nphi, Nphi);
C = sparse(Nphi, Npatches);

% indexing for phi terms (in mV):
% phi_m,i:      (i-1)*(M+3)+1,          i = 1:N
% phi_R,i:      (i-1)*(M+3)+2,          i = 1:N-1
% phi_L,i:      (i-1)*(M+3),            i = 2:N
% phi_c,i,j:    (i-1)*(M+3)+j+2,        i = 1:N-1, j = 1:M
%
% indexing for I terms (in uA):
% I_m,i:        (i-1)*(2*M+1)+1,        i = 1:N
% I_R,i,j:      (i-1)*(2*M+1)+j+1,      i = 1:N-1, j = 1:M
% I_L,i,j:      (i-2)*(2*M+1)+M+j+1,    i = 2:N,   j = 1:M

for i = 1:N
    
    % c (column) terms:
    c_phi_m = (i-1)*(M+3)+1;         % phi_m,i
    c_phi_R = (i-1)*(M+3)+2;         % phi_R,i
    c_phi_Rm1 = (i-1-1)*(M+3)+2;     % phi_R,i-1
    c_phi_L = (i-1)*(M+3);           % phi_L,i
    c_phi_Lp1 = (i-1+1)*(M+3);       % phi_L,i+1
    
    c_I_m = (i-1)*(2*M+1)+1;         % I_m,i
    
    % EQN. 1, rows 1:N
    r = i;  % row
    if i == 1
        P(r,c_phi_m) = (Cm/dt + gmyo);
        P(r,c_phi_R) = -gmyo;
        Q(r,c_phi_m) = Cm/dt;
        C(r,c_I_m) = -1;
    elseif i == N
        P(r,c_phi_m) = (Cm/dt + gmyo);
        P(r,c_phi_L) = -gmyo;
        Q(r,c_phi_m) = Cm/dt;
        C(r,c_I_m) = -1;
    else
        P(r,c_phi_m) = (Cm/dt + 2*gmyo);
        P(r,c_phi_L) = -gmyo;
        P(r,c_phi_R) = -gmyo;
        Q(r,c_phi_m) = Cm/dt;
        C(r,c_I_m) = -1;
    end
    
    % EQN. 2, rows N+1:2*N-1
    if i > 1
        r = i + N - 1;
        P(r,c_phi_m) = -gmyo;
        P(r,c_phi_L) = (gmyo + ggap(i-1) + sum(CL(i,:)/dt));
        P(r,c_phi_Rm1) = -ggap(i-1);
        Q(r,c_phi_L) = sum(CL(i,:)/dt);
        for j = 1:M
            c_phi_cm1 = (i-1-1)*(M+3)+j+2;       % phi_c,i-1,j
            c_I_L = (i-2)*(2*M+1)+M+j+1;         % I_L,i,j
            
            P(r,c_phi_cm1) = -CL(i,j)/dt;
            Q(r,c_phi_cm1) = -CL(i,j)/dt;
            C(r,c_I_L) = -1;
        end
    end
    
    % EQN. 3, rows 2*N:3*N-2
    if i < N
        r = i + 2*N - 1;
        P(r,c_phi_m) = -gmyo;
        P(r,c_phi_Lp1) = -ggap(i);
        P(r,c_phi_R) = (gmyo + ggap(i) + sum(CR(i,:)/dt));
        Q(r,c_phi_R) = sum(CR(i,:)/dt);
        for j = 1:M
            c_phi_c = (i-1)*(M+3)+j+2;     % phi_c,i,j
            c_I_R = (i-1)*(2*M+1)+j+1;     % I_R,i,j
            
            P(r,c_phi_c) = -CR(i,j)/dt;
            Q(r,c_phi_c) = -CR(i,j)/dt;
            C(r,c_I_R) = -1;
        end
    end
    
    % EQN. 4, rows 3*N-1:3*N-2+(N-1)*M
    if i < N
        for j = 1:M
            r = j + (i-1)*M + 3*N - 2;  %
            knoteqj = setdiff(1:M,j);
            c_phi_c = (i-1)*(M+3)+j+2;        % phi_c,i,j
            c_I_Lp1 = (i-2+1)*(2*M+1)+M+j+1;  % I_L,i+1,j
            c_I_R = (i-1)*(2*M+1)+j+1;        % I_R,i,j
            
            P(r,c_phi_Lp1) = CL(i+1,j)/dt;
            P(r,c_phi_R) = CR(i,j)/dt;
            P(r,c_phi_c) = (-ns*Gb(i,j) + -ns*sum(Gc(j,knoteqj,i)) - CR(i,j)/dt - CL(i+1,j)/dt);
            Q(r,c_phi_Lp1) = CL(i+1,j)/dt;
            Q(r,c_phi_R) = CR(i,j)/dt;
            Q(r,c_phi_c) = (-CR(i,j)/dt - CL(i+1,j)/dt);
            C(r,c_I_Lp1) = -1;
            C(r,c_I_R) = -1;
            for k = knoteqj
                c_phi_ck = (i-1)*(M+3)+k+2;    % phi_c,i,k
                P(r,c_phi_ck) = ns*Gc(j,k,i);
            end
            
        end
    end
    
    
    
end
