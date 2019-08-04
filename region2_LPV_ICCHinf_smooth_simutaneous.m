%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code is for the DOF mixed ICC/Hinf LPV controller design for Active Magnetic
%%% Bearing system. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
rolmip('clearvar')   % clear all variables in workspace; see fucntion 'rolmip' for details 
clear , clc
warning('off','YALMIP:strict');
% global switchflag;
% switchflag = 1;
% global A_LPV1;
%% formulate orginal system 
%%%%%% form state-space model of AMB 

parabound = [315, 1100]; %%% 
range = parabound(2)-parabound(1);

parabound1 = [315,720];
range1 = parabound1(2) - parabound1(1);

parabound2 = [700,1100];
range2 = parabound2(2) - parabound2(1);

kappa = 100;  %%% varying rate bound of scheduling parameter
pararatebound = [-kappa, kappa]; 

%%%%% Formulate system matrix for entire parameter region %%%%%
m =  14.46; %%% unit kg
Ja = 0.0136;  %%% kg*m^2
Jr = 0.333;  %%%kg*m^2
leng = 0.13; %%% unit: m 
m1 = Jr/(leng^2); 
N = 400 ;   %%% turns
k = 4.68e8; 
R = 14.7; %%%Omega
mu0 = 4*pi*10^-7; %%% H/m; permeability of free space.
Phi_0 = 2.09e-4; %%% Wb
G_0 = 0.55e-3; %%% meter
h = 40e-3; %%% meter
A = 1531.79e-6 ;%%% meter*meter; the area under one electromagnet pole (each electromagnet has two),
c1 =  2*k*Phi_0*(1+2*G_0/(pi*h));
c2 =  2*k*((Phi_0)^2)/(pi*h);
d1 = 2*R*G_0/(mu0*A*N);
d2 = 2*R*Phi_0/(mu0*A*N);

Ap1 = [ 0,    0,    1,    0,     0,     0; 
       0,    0,    0,    1,     0,     0;
       -4*c2/m1, 0 ,0,    0,     2*c1/m1, 0;
       0, -4*c2/m1,  0,   0,     0,     2*c1/m1;
       2*d2/N,0,   0,    0,     -d1/N,   0;
       0,   2*d2/N, 0,   0,     0,      -d1/N ];
   
Ap2 = [0, 0,  0, 0 ,0,0;
       0, 0, 0, 0, 0, 0;
       0,  0, 0, -Ja/Jr, 0 ,0;
       0 , 0, Ja/Jr, 0 ,0 ,0;
       0, 0,  0, 0 ,0,0;
       0, 0, 0, 0, 0, 0];
   
Ap{1} = {0,Ap1};
Ap{2} = {1,Ap2};
%%% formulate A_LPV1 in parabound1
Ap_LPV = rolmipvar(Ap,'Ap_LPV',parabound);

Ap_LPV1 = rolmipvar(Ap,'Ap_LPV1',parabound1);
Ap_LPV2 = rolmipvar(Ap,'Ap_LPV2',parabound2);

Bu = (1/N)*[zeros(4,2);eye(2)];
Bw = zeros(6,2);
% B1 = [zeros(2,2);eye(2);zeros(2,2)]/m1; % All model share same Binf
%%% assume all disturbance has same influence on system

% Czinf = [eye(2),zeros(2,4)];
% Cz2 = [eye(2),zeros(2,4)];
Cp = [eye(2),zeros(2,4)];

%%%% D matrix  %%%%

 
Dzinfw = zeros(2,2);
Dyw = eye(2);
Dyu = zeros(2,2);

%%%%%% weighting functions %%%%%%

Wz_1 = tf([1e1 8e1],[1 1e-3]);
Wz = [Wz_1]*eye(2);
Wz_sys = ss(Wz,'minimal');

Wu_1 = tf([0.01 1],[1 1e5]);
Wu = [Wu_1]*eye(2);
Wu_sys = ss(Wu,'minimal');


[m_Az,n_Az] = size(Wz_sys.A);
[m_Au,n_Au] = size(Wu_sys.A);
[m_Cz,n_Cz] = size(Wz_sys.C);
[m_Cu,n_Cu] = size(Wu_sys.C);

% %%%%%% augmented system matrix %%%%%
[m_Ap,n_Ap] = size(Ap_LPV);
[m_Cz2,n_Cz2] = size(Cp);
[m_Dzinfw,n_Dzinfw] = size(Dzinfw);
%%%% 
A_LPV = [Ap_LPV,zeros(m_Ap,n_Az),zeros(m_Ap,n_Au);...
         Wz_sys.B*Cp,Wz_sys.A,zeros(m_Az,n_Au);...
         zeros(m_Au,n_Ap),zeros(m_Au,n_Az),Wu_sys.A];
     
B1 = [Bw;zeros(2,2);zeros(m_Au,2)];
B2 = [Bu;zeros(2,2);Wu_sys.B];

A_LPV1 = [Ap_LPV1,zeros(m_Ap,n_Az),zeros(m_Ap,n_Au);...
         Wz_sys.B*Cp,Wz_sys.A,zeros(m_Az,n_Au);...
         zeros(m_Au,n_Ap),zeros(m_Au,n_Az),Wu_sys.A];

A_LPV2 = [Ap_LPV2,zeros(m_Ap,n_Az),zeros(m_Ap,n_Au);...
         Wz_sys.B*Cp,Wz_sys.A,zeros(m_Az,n_Au);...
         zeros(m_Au,n_Ap),zeros(m_Au,n_Az),Wu_sys.A];
%%%

C1 = [Wz_sys.D*Cp, Wz_sys.C, zeros(m_Cz,n_Cu);...
      zeros(m_Cu,n_Ap),zeros(m_Cu,n_Cz),Wu_sys.C];
D11 = [zeros(2,2); zeros(2,2)];
D12 = [zeros(2,2);Wu_sys.D];

%  C2 = [Cp,zeros(2,2),zeros(m_Cz,2)];
C2 = [Wz_sys.D*Cp, Wz_sys.C, zeros(m_Cz,n_Cu)];

E1 = zeros(2,2);
E2 = zeros(2,2);

Cy = [Cp,zeros(m_Cz2,n_Az),zeros(m_Cz2,n_Au)];
Dy = [0.001*Dyw];

%%%%%%%

[m_A,n_A] = size(A_LPV);
[m_B2,n_B2] = size(B2);
[m_B1,n_B1] = size(B1);
[m_C1,n_C1] = size(C1);
[m_C2,n_C2] = size(C2); 
[m_Cy,n_Cy] = size(Cy); 
[m_D11,n_D11] = size(D11);
[m_D12,n_D12] = size(D12);
[m_Dy,n_Dy] = size(Dy);
[m_E2,n_E2] = size(E2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
R = rolmipvar(m_A,n_A,'R','symmetric',0,0); %%% set X = constant 
% X_inf = rolmipvar(m_A,n_A,'X_inf','symmetric',{[0],[1]},parabound);
dotR  = rolmipdiff(R,'dotR2',parabound,pararatebound);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Y_inf = rolmipvar(m_A,n_A,'Y_inf','symmetric',0,0); %%% set Y = constant 
S1 = rolmipvar(m_A,n_A,'S1','symmetric',{0,1},parabound1);
dotS1  =  rolmipdiff(S1,'dotS1',parabound1,pararatebound); 

S2 = rolmipvar(m_A,n_A,'S2','symmetric',{0,1},parabound2);
dotS2  =  rolmipdiff(S2,'dotS2',parabound2,pararatebound); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = rolmipvar(m_C2,m_C2,'W1','symmetric',0,0);

% W = rolmipvar(m_C2,m_C2,'W2','symmetric',0,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hat_Ak1 = rolmipvar(m_A,n_A,'hat_Ak1','full',{0,1},parabound1);
hat_Bk1 = rolmipvar(m_A,m_Cy,'hat_Bk1','full',{0,1},parabound1);
hat_Ck1 = rolmipvar(n_B2,n_A,'hat_Ck1','full',{0,1},parabound1);

hat_Ak2 = rolmipvar(m_A,n_A,'hat_Ak2','full',{0,1},parabound2);
hat_Bk2 = rolmipvar(m_A,m_Cy,'hat_Bk2','full',{0,1},parabound2);
hat_Ck2 = rolmipvar(n_B2,n_A,'hat_Ck2','full',{0,1},parabound2);

%%% 
mu = rolmipvar(36,'mu','scalar'); 
%  gamma = sdpvar(1);
TRW = sdpvar(1);

Ubar = 10e7;

beta = 3; 
% beta = 2;

%%%%%%%%%%%%%    LMIs  %%%%%%%%%%%%%%%%%%%%%%

Lmi11{11} = -dotR + (A_LPV1*R + B2*hat_Ck1) + (A_LPV1*R + B2*hat_Ck1)';
Lmi21{11} = hat_Ak1 + A_LPV1' ;
Lmi31{11} = B1';
Lmi41{11} = C1*R + D12*hat_Ck1;
Lmi22{11} = (S1*A_LPV1 + hat_Bk1*Cy) + (S1*A_LPV1 + hat_Bk1*Cy)' + dotS1;
Lmi32{11} = (S1*B1 + hat_Bk1*Dy)'; 
Lmi42{11} = C1;
Lmi33{11} = -eye(n_D11);
Lmi43{11} = D11;
Lmi44{11} = -mu*eye(m_D11);


T11 = [Lmi11{11},(Lmi21{11})',(Lmi31{11})',(Lmi41{11})';
     Lmi21{11},Lmi22{11},(Lmi32{11})',(Lmi42{11})';
     Lmi31{11},Lmi32{11},Lmi33{11},(Lmi43{11})';
     Lmi41{11},Lmi42{11},Lmi43{11},Lmi44{11}];

T12 = [R,eye(m_A);eye(m_A),S1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lmi11{21} = -dotR + (A_LPV2*R + B2*hat_Ck2) + (A_LPV2*R + B2*hat_Ck2)';
Lmi21{21} = hat_Ak2 + A_LPV2' ;
Lmi31{21} = B1';
Lmi41{21} = C1*R + D12*hat_Ck2;
Lmi22{21} = (S2*A_LPV2 + hat_Bk2*Cy) + (S2*A_LPV2 + hat_Bk2*Cy)' + dotS2;
Lmi32{21} = (S2*B1 + hat_Bk2*Dy)'; 
Lmi42{21} = C1;
Lmi33{21} = -eye(n_D11);
Lmi43{21} = D11;
Lmi44{21} = -mu*eye(m_D11);


T21 = [Lmi11{21},(Lmi21{21})',(Lmi31{21})',(Lmi41{21})';
     Lmi21{21},Lmi22{21},(Lmi32{21})',(Lmi42{21})';
     Lmi31{21},Lmi32{21},Lmi33{21},(Lmi43{21})';
     Lmi41{21},Lmi42{21},Lmi43{21},Lmi44{21}];

T22 = [R,eye(m_A);eye(m_A),S2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lmi11{12} = -dotR + (A_LPV1*R + B2*hat_Ck1) + (A_LPV1*R + B2*hat_Ck1)';
Lmi21{12} = hat_Ak1 + A_LPV1' ;
Lmi31{12} = B1';
Lmi22{12} = (S1*A_LPV1 + hat_Bk1*Cy) + (S1*A_LPV1 + hat_Bk1*Cy)' + dotS1;
Lmi32{12} = (S1*B1 + hat_Bk1*Dy)';
Lmi33{12} = -eye(n_D11);

T13 = [Lmi11{12},(Lmi21{12})',(Lmi31{12})';
     Lmi21{12},Lmi22{12},(Lmi32{12})';
     Lmi31{12},Lmi32{12},Lmi33{12}];

Lmi11{13} = W;
Lmi21{13} = (C2*R)';
Lmi31{13} = C2';
Lmi22{13} = R ;
Lmi32{13} = eye(n_A);
Lmi33{13} = S1;

T14 = [Lmi11{13},(Lmi21{13})',(Lmi31{13})';
     Lmi21{13},Lmi22{13},(Lmi32{13})';
     Lmi31{13},Lmi32{13},Lmi33{13}]; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lmi11{22} = -dotR + (A_LPV2*R + B2*hat_Ck2) + (A_LPV2*R + B2*hat_Ck2)';
Lmi21{22} = hat_Ak2 + A_LPV2' ;
Lmi31{22} = B1';
Lmi22{22} = (S2*A_LPV2 + hat_Bk2*Cy) + (S2*A_LPV2 + hat_Bk2*Cy)' + dotS2;
Lmi32{22} = (S2*B1 + hat_Bk2*Dy)';
Lmi33{22} = -eye(n_D11);

T23 = [Lmi11{22},(Lmi21{22})',(Lmi31{22})';
     Lmi21{22},Lmi22{22},(Lmi32{22})';
     Lmi31{22},Lmi32{22},Lmi33{22}];

Lmi11{23} = W;
Lmi21{23} = (C2*R)';
Lmi31{23} = C2';
Lmi22{23} = R ;
Lmi32{23} = eye(n_A);
Lmi33{23} = S2;

T24 = [Lmi11{23},(Lmi21{23})',(Lmi31{23})';
     Lmi21{23},Lmi22{23},(Lmi32{23})';
     Lmi31{23},Lmi32{23},Lmi33{23}]; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi{1} = [1,0]; 
Phi{2} = [0,1];
 
Lmi11{14} = Ubar ;
Lmi21{14} = (Phi{1}*hat_Ck1)';
Lmi31{14} = zeros(1,n_A)';
Lmi22{14} = R ;
Lmi32{14} = eye(n_A);
Lmi33{14} = S1;

T15 = [Lmi11{14},(Lmi21{14})',(Lmi31{14})';
     Lmi21{14},Lmi22{14},(Lmi32{14})';
     Lmi31{14},Lmi32{14},Lmi33{14}]; 

Lmi11{15} = Ubar ;
Lmi21{15} = (Phi{2}*hat_Ck1)';
Lmi31{15} = zeros(1,n_A)';
Lmi22{15} = R ;
Lmi32{15} = eye(n_A);
Lmi33{15} = S1;

T16 = [Lmi11{15},(Lmi21{15})',(Lmi31{15})';
     Lmi21{15},Lmi22{15},(Lmi32{15})';
     Lmi31{15},Lmi32{15},Lmi33{15}]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lmi11{24} = Ubar ;
Lmi21{24} = (Phi{1}*hat_Ck2)';
Lmi31{24} = zeros(1,n_A)';
Lmi22{24} = R ;
Lmi32{24} = eye(n_A);
Lmi33{24} = S2;

T25 = [Lmi11{24},(Lmi21{24})',(Lmi31{24})';
     Lmi21{24},Lmi22{24},(Lmi32{24})';
     Lmi31{24},Lmi32{24},Lmi33{24}]; 

Lmi11{25} = Ubar ;
Lmi21{25} = (Phi{2}*hat_Ck2)';
Lmi31{25} = zeros(1,n_A)';
Lmi22{25} = R ;
Lmi32{25} = eye(n_A);
Lmi33{25} = S2;

T26 = [Lmi11{25},(Lmi21{25})',(Lmi31{25})';
     Lmi21{25},Lmi22{25},(Lmi32{25})';
     Lmi31{25},Lmi32{25},Lmi33{25}]; 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S2_S12 = (evalpar(S2,{[1-(parabound1(2)-parabound2(1))/range2,(parabound1(2)-parabound2(1))/range2]}));
%%%%%
S2_S21 = (evalpar(S2,{[1,0]}));
%%%%%
S1_S12 = evalpar(S1,{[0,1]});
%%%%%
S1_S21 = evalpar(S1,{[1-(parabound2(1)-parabound1(1))/range1,(parabound2(1)-parabound1(1))/range1]});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hat_Ak2_S12 = (evalpar(hat_Ak2,{[1-(parabound1(2)-parabound2(1))/range2,(parabound1(2)-parabound2(1))/range2]}));
%%%%%
hat_Ak2_S21 = (evalpar(hat_Ak2,{[1,0]}));
%%%%%
hat_Ak1_S12 = evalpar(hat_Ak1,{[0,1]});
%%%%%
hat_Ak1_S21 = evalpar(hat_Ak1,{[1-(parabound2(1)-parabound1(1))/range1,(parabound2(1)-parabound1(1))/range1]});


hat_Bk2_S12 = (evalpar(hat_Bk2,{[1-(parabound1(2)-parabound2(1))/range2,(parabound1(2)-parabound2(1))/range2]}));
%%%%%
hat_Bk2_S21 = (evalpar(hat_Bk2,{[1,0]}));
%%%%%
hat_Bk1_S12 = evalpar(hat_Bk1,{[0,1]});
%%%%%
hat_Bk1_S21 = evalpar(hat_Bk1,{[1-(parabound2(1)-parabound1(1))/range1,(parabound2(1)-parabound1(1))/range1]});

hat_Ck2_S12 = (evalpar(hat_Ck2,{[1-(parabound1(2)-parabound2(1))/range2,(parabound1(2)-parabound2(1))/range2]}));
%%%%%
hat_Ck2_S21 = (evalpar(hat_Ck2,{[1,0]}));
%%%%%
hat_Ck1_S12 = evalpar(hat_Ck1,{[0,1]});
%%%%%
hat_Ck1_S21 = evalpar(hat_Ck1,{[1-(parabound2(1)-parabound1(1))/range1,(parabound2(1)-parabound1(1))/range1]});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% LMI1 = [T11 <= 0; T12 >= 0; T13 <= 0 ; T14 >= 0 ; T15 >= 0 ; T16 >= 0; ];
% LMI2 = [T21 <= 0; T22 >= 0; T23 <= 0 ; T24 >= 0 ; T25 >= 0 ; T26 >= 0; ];

LMI1 = [T11 <= 0; T12 >= 0; T13 <= 0 ; T14 >= 0 ;];
LMI2 = [T21 <= 0; T22 >= 0; T23 <= 0 ; T24 >= 0 ;];


LMI3 = [S1_S12 >= S2_S12; S2_S21 >= S1_S21  ];
LMI4 = [trace(W) <= TRW; TRW >= 0 ]; 

LMIs = LMI1 + LMI2 + LMI3 + LMI4; 


norm_swit_S12 = norm(hat_Ak1_S12 - hat_Ak2_S12) + norm(hat_Bk1_S12 - hat_Bk2_S12)...
    + norm(hat_Ck1_S12 - hat_Ck2_S12) + norm(S1_S12 - S2_S12);

norm_swit_S21 = norm(hat_Ak1_S21 - hat_Ak2_S21) + norm(hat_Bk1_S21 - hat_Bk2_S21)...
    + norm(hat_Ck1_S21 - hat_Ck2_S21) + norm(S1_S21 - S2_S21);


% trW = vertices(trace(W));
% obj = [trW{1}];
% obj = [gamma];

% obj = (10^(beta)*TRW +(norm_swit_S12 + norm_swit_S21));
% obj = (beta*TRW);
obj = [(1e8)*TRW];
% obj = [];

%%%%%  

% Sol = optimize(LMIs,obj,sdpsettings('solver','sedumi','verbose',0));
% Sol = solvesdp(LMIs,obj,sdpsettings('solver','sedumi','debug',1));
Sol = solvesdp(LMIs,obj,sdpsettings('solver','sedumi','verbose',0));

Info.gamma = sqrt(double(mu));
Info.TraceW = double(TRW);
Info.smoothness = double(norm_swit_S12 + norm_swit_S21);


if Sol.problem == 0
    Info.Solution = 'Successfully solved';
    elseif Sol.problem == 4
    Info.Solution = 'Numerical problems';
    elseif Sol.problem == 1
    Info.Solution = 'Infeasible';
    elseif Sol.problem == -4
    Info.Solution = 'Solver not applicable';
end

Info


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
grid_step = 5 ; % define gridding step  
para = parabound(1):grid_step:parabound(2);
para1 = parabound1(1):grid_step:parabound1(2);
para2 = parabound2(1):grid_step:parabound2(2);

% %%% formulate LPV model and controller %%% 

sysB = [B1,B2];
sysC = [C1;C2;Cy];
sysD = [D11,D12;E1,E2;Dy,Dyu];

% %%%%% extract and construct LPV controller %%%%%
% 
R0 = double(evalpar(R,{[1,0]}));
invM_T0 = eye(m_A);

Dk = zeros(2,2);
 
for grid_i = 1:length(para)
    
  sysA{grid_i} = double(evalpar(A_LPV,{[1-(para(grid_i)-parabound(1))/range,(para(grid_i)-parabound(1))/range]}));
  sys(:,:,grid_i) = ss(sysA{grid_i},sysB,sysC,sysD);    
  EIG{grid_i} = eig(sysA{grid_i});

%   CTRB_rank{i} = rank(ctrb(sysA{i},B2));
%   OBSV_rank{i} = rank(obsv(sysA{i},Cy));
%   
%   [K{i},CL{i},GAM{i},INFO{i}] = hinfsyn(sys1(:,:,i),2,2,'method','lmi');
    
%   open_sysA{i} = double(evalpar(Ap_LPV,{[1-(para(i)-parabound(1))/range,(para(i)-parabound(1))/range]}));
%   CTRB_rank_open{i} = rank(ctrb(open_sysA{i},Bu));
%   
%   Ctrb_unstable_mode{i} = rank([EIG{i}(3,1)*eye(10)-sysA{i},B2]);

end
sys.SamplingGrid = struct('para',para);


sysA1_vert1 = double(evalpar(A_LPV1,{[1,0]}));
sysA1_vert2 = double(evalpar(A_LPV1,{[0,1]}));

S1_vert1 = double(evalpar(S1,{[1,0]}));
S1_vert2 = double(evalpar(S1,{[0,1]}));

hatAk1_vert1 = double(evalpar(hat_Ak1,{[1,0]}));
hatAk1_vert2 = double(evalpar(hat_Ak1,{[0,1]}));

hatBk1_vert1 = double(evalpar(hat_Bk1,{[1,0]}));
hatBk1_vert2 = double(evalpar(hat_Bk1,{[0,1]}));

hatCk1_vert1 = double(evalpar(hat_Ck1,{[1,0]}));
hatCk1_vert2 = double(evalpar(hat_Ck1,{[0,1]}));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sysA2_vert1 = double(evalpar(A_LPV2,{[1,0]}));
sysA2_vert2 = double(evalpar(A_LPV2,{[0,1]}));

S2_vert1 = double(evalpar(S2,{[1,0]}));
S2_vert2 = double(evalpar(S2,{[0,1]}));

hatAk2_vert1 = double(evalpar(hat_Ak2,{[1,0]}));
hatAk2_vert2 = double(evalpar(hat_Ak2,{[0,1]}));

hatBk2_vert1 = double(evalpar(hat_Bk2,{[1,0]}));
hatBk2_vert2 = double(evalpar(hat_Bk2,{[0,1]}));

hatCk2_vert1 = double(evalpar(hat_Ck2,{[1,0]}));
hatCk2_vert2 = double(evalpar(hat_Ck2,{[0,1]}));


% for grid_i1 = 1:length(para1)
% 
%   sysA1{grid_i1} = double(evalpar(A_LPV1,{[1-(para1(grid_i1)-parabound1(1))/range1,(para1(grid_i1)-parabound1(1))/range1]}));
% 
%   S1_grid{grid_i1} = double(evalpar(S1,{[1-(para1(grid_i1)-parabound1(1))/range1,(para1(grid_i1)-parabound1(1))/range1]}));
%   N1_grid{grid_i1} = (eye(m_A) - S1_grid{grid_i1}*R0);
%   U1_grid{grid_i1} = -S1_grid{grid_i1}*(inv(N1_grid{grid_i1}))';
%   
%   hatAk1_grid{grid_i1} = double(evalpar(hat_Ak1,{[1-(para1(grid_i1)-parabound1(1))/range1,(para1(grid_i1)-parabound1(1))/range1]}));
%   hatBk1_grid{grid_i1} = double(evalpar(hat_Bk1,{[1-(para1(grid_i1)-parabound1(1))/range1,(para1(grid_i1)-parabound1(1))/range1]}));
%   hatCk1_grid{grid_i1} = double(evalpar(hat_Ck1,{[1-(para1(grid_i1)-parabound1(1))/range1,(para1(grid_i1)-parabound1(1))/range1]}));
%   
%   Ak1{grid_i1} = inv(N1_grid{grid_i1})*(hatAk1_grid{grid_i1} - S1_grid{grid_i1}*sysA1{grid_i1}*R0 - hatBk1_grid{grid_i1}*Cy*R0 - S1_grid{grid_i1}*B2*hatCk1_grid{grid_i1})*invM_T0;
%   
%   Bk1{grid_i1} = inv(N1_grid{grid_i1})*(hatBk1_grid{grid_i1});
%   
%   Ck1{grid_i1} = (hatCk1_grid{grid_i1})*invM_T0;
%   
%   
%   ICC_cost1_U1{grid_i1} = [1,0]*Ck1{grid_i1}*U1_grid{grid_i1}*(Ck1{grid_i1})'*[1;0];
%   
%   ICC_cost1_U2{grid_i1} = [0,1]*Ck1{grid_i1}*U1_grid{grid_i1}*(Ck1{grid_i1})'*[0;1];
%   
%   controller1(:,:,grid_i1) = ss(Ak1{grid_i1},Bk1{grid_i1},Ck1{grid_i1},Dk);
%   
%   A_cl1{grid_i1} = [sysA1{grid_i1},B2*Ck1{grid_i1};Bk1{grid_i1}*Cy,Ak1{grid_i1}];
%   B_cl1{grid_i1} = [B1;Bk1{grid_i1}*Dy];
%   C_inf_cl1{grid_i1} = [C1,D12*Ck1{grid_i1}];
%   C_2_cl1{grid_i1} = [C2,E2*Ck1{grid_i1}];
% %   TRCPC{i} = trace(C2*R0*C2');
%   EIG_cl1{grid_i1} = eig(A_cl1{grid_i1});
%   
% end
% 
% % sys_open.SamplingGrid = struct('para',para);
% controller1.SamplingGrid = struct('para1',para1);
% 
% 
% for grid_i2 = 1:length(para2)
% 
%   sysA2{grid_i2} = double(evalpar(A_LPV2,{[1-(para2(grid_i2)-parabound2(1))/range2,(para2(grid_i2)-parabound2(1))/range2]}));
% 
%   S2_grid{grid_i2} = double(evalpar(S2,{[1-(para2(grid_i2)-parabound2(1))/range2,(para2(grid_i2)-parabound2(1))/range2]}));
%   N2_grid{grid_i2} = (eye(m_A) - S2_grid{grid_i2}*R0);
%   U2_grid{grid_i2} = -S2_grid{grid_i2}*(inv(N2_grid{grid_i2}))';
%   
%   hatAk2_grid{grid_i2} = double(evalpar(hat_Ak2,{[1-(para2(grid_i2)-parabound2(1))/range2,(para2(grid_i2)-parabound2(1))/range2]}));
%   hatBk2_grid{grid_i2} = double(evalpar(hat_Bk2,{[1-(para2(grid_i2)-parabound2(1))/range2,(para2(grid_i2)-parabound2(1))/range2]}));
%   hatCk2_grid{grid_i2} = double(evalpar(hat_Ck2,{[1-(para2(grid_i2)-parabound2(1))/range2,(para2(grid_i2)-parabound2(1))/range2]}));
%   
%   Ak2{grid_i2} = inv(N2_grid{grid_i2})*(hatAk2_grid{grid_i2} - S2_grid{grid_i2}*sysA2{grid_i2}*R0 - hatBk2_grid{grid_i2}*Cy*R0 - S2_grid{grid_i2}*B2*hatCk2_grid{grid_i2})*invM_T0;
%   
%   Bk2{grid_i2} = inv(N2_grid{grid_i2})*(hatBk2_grid{grid_i2});
%   
%   Ck2{grid_i2} = (hatCk2_grid{grid_i2})*invM_T0;
%   
%   
%   ICC_cost2_U1{grid_i2} = [1,0]*Ck2{grid_i2}*U2_grid{grid_i2}*(Ck2{grid_i2})'*[1;0];
%   
%   ICC_cost2_U2{grid_i2} = [0,1]*Ck2{grid_i2}*U2_grid{grid_i2}*(Ck2{grid_i2})'*[0;1];
%   
%   controller2(:,:,grid_i2) = ss(Ak2{grid_i2},Bk2{grid_i2},Ck2{grid_i2},Dk);
%   
%   A_cl2{grid_i2} = [sysA2{grid_i2},B2*Ck2{grid_i2};Bk2{grid_i2}*Cy,Ak2{grid_i2}];
%   B_cl2{grid_i2} = [B1;Bk2{grid_i2}*Dy];
%   C_inf_cl2{grid_i2} = [C1,D12*Ck2{grid_i2}];
%   C_2_cl2{grid_i2} = [C2,E2*Ck2{grid_i2}];
% %   TRCPC{i} = trace(C2*R0*C2');
%   EIG_cl2{grid_i2} = eig(A_cl2{grid_i2});
%     
% end
% 
% % sys_open.SamplingGrid = struct('para',para);
% controller2.SamplingGrid = struct('para2',para2);

control_state = zeros(10,1);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% grid_num = 1:1:length(para); 
% 
% for num = grid_num(1):grid_num(end)
%  
% EIG_real{num} = real(EIG{num});
% EIG_img{num} = imag(EIG{num});
% 
% EIGcl_real{num} = real(EIG_cl{num});
% EIGcl_img{num} = imag(EIG_cl{num});
% 
% end
% 
% 
% figure(1)
% for thetaset = grid_num(1):grid_num(end)
% ax = plot(EIG_real{thetaset},EIG_img{thetaset},'^',...
%     'MarkerSize',7,...
%      'MarkerEdgeColor',(1/length(grid_num))*(thetaset-grid_num(1))*[0,0,1],...
%     'MarkerFaceColor',(1/length(grid_num))*(thetaset-grid_num(1))*[0,0,1]);
% hold on
% 
% plot(EIGcl_real{thetaset},EIGcl_img{thetaset},'o',...
%     'MarkerSize',7,...
%      'MarkerEdgeColor',(1/length(grid_num))*(thetaset-grid_num(1))*[1,0,0],...
%     'MarkerFaceColor',(1/length(grid_num))*(thetaset-grid_num(1))*[1,0,0]);
% hold on
% 
% end
% 
% 
% colorbar
% grid on
% % axis([-1000,100,-200,200])
% xlabel('Real')
% ylabel('Imag')
% title('open-loop closed-loop root locus')
