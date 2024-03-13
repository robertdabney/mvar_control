% ========================================================================
% PLANT:
% ========================================================================
p = @(varargin)1/sum(1./[varargin{:}]);
syms s;
% s = tf([1,0],1);

L1=.001;
ZL1=s*L1;
L2=.05;
ZL2=s*L2;
L3=.01;
ZL3=s*L3;
C1=.01;
ZC1=1/(s*C1);
C2=.005;
ZC2=1/(s*C2);
C3=.001;
ZC3=1/(s*C3);
R1=100;
ZR1=R1;
R2=100;
ZR2=R2;

%====G11====
Z1=p(ZC1+ZL1,ZR1);
Z2=ZR2+ZL2;
Z3=p(ZL3+ZC3,ZC2);
Ztot = p(Z2,Z3)/Z2;
G11 = Ztot/(p(Z2,Z3)+Z1);
[Gnumsym,Gdensym] = numden(G11);
Gnum=sym2poly(Gnumsym);
Gden=sym2poly(Gdensym);
G11=tf(Gnum,Gden);

%====G12====
Z1 = ZC1+ZL1;
Z2 = p(p(ZR2+ZL2,ZC2),ZL3+ZC3);
G12 = Z1/((ZR1+Z1)*(p(Z1,ZR1)+Z2));
[Gnumsym,Gdensym] = numden(G12);
Gnum=sym2poly(Gnumsym);
Gden=sym2poly(Gdensym);
G12=tf(Gnum,Gden);

%====G21====
Z1 = p((ZC1+ZL1),ZR1);
Z2 = ZR2+ZL2;
Z3 = p(ZC2,(ZL3+ZC3));
G21 = -1/(p(Z1,Z3)+Z2);
[Gnumsym,Gdensym] = numden(G21);
Gnum=sym2poly(Gnumsym);
Gden=sym2poly(Gdensym);
G21=tf(Gnum,Gden);

%====G22====
Z1 = ZC1+ZL1;
% Z2 and Z3 are the same as G21
Ztot = p(ZR1,p(Z1,Z3))+Z2;
G22 = -(p(Z1,Z3)/(ZR1+p(Z1,Z3)))/Ztot;
[Gnumsym,Gdensym] = numden(G22);
Gnum=sym2poly(Gnumsym);
Gden=sym2poly(Gdensym);
G22=tf(Gnum,Gden);


G_tf =[G11,G21;G12,G22];
G_ss = ss(G_tf);

% Minimized Realization
G_min_ss = minreal(G_ss);

% Smith Mc-Millan Form
[SM,U,V] = smform(G_min_ss);
% Coprime Factorization
[fact,Mr,Nr] = rncf(G_tf);

P22_poles = pole(G_min_ss);
P22_tzeros = tzero(G_min_ss);
G11_zeros = zero(G_min_ss(1,1));
G11_poles = pole(G_min_ss(1,1));
G12_zeros = zero(G_min_ss(1,2));
G12_poles = pole(G_min_ss(1,2));
G21_zeros = zero(G_min_ss(2,1));
G21_poles = pole(G_min_ss(2,1));
G22_zeros = zero(G_min_ss(2,2));
G22_poles = pole(G_min_ss(2,2));

num_states = size(G_ss.A,1);
num_u_in = 2;
num_v_out = 2;
num_w_in = 4;
num_z_out = 3;

epsilon1 = 1;
epsilon2 = 1;

B1 = [zeros(num_states,num_w_in-1),G_ss.B(:,1)];
C1 = [G_ss.C(1,:);zeros(num_z_out-1,num_states)];
D11 = zeros(num_z_out,num_w_in);D11(1,1)=-1;
D12 = zeros(num_z_out,num_u_in);D12(2,1)=epsilon1;D12(3,2)=epsilon2;
D21 = [-1,1,0,0;0,0,1,0];
B = [B1,G_ss.B];
C = [C1;G_ss.C];
D = [D11,D12;D21,G_ss.D];
P = ss(G_ss.A,B,C,D);
B2 = G_ss.B;
C2 = G_ss.C;
D22 = G_ss.D;
P11 = ss(G_ss.A,B1,C1,D11);
P12 = ss(G_ss.A,B2,C1,D12);
P21 = ss(G_ss.A,B1,C2,D21);
P22 = ss(G_ss.A,B2,C2,D22);

Q = eye(2);
K = q_control(G_ss,Q);


closed_loop_robustness = feedback(G_ss*K,eye(2));
figure(1)
% SISO nyquist
nyquist(G_ss*K);
w=logspace(-2,2,100);
% MIMO Nyquist
gamma=mvar_nyquist(G_ss*K,w);
figure(2); 
plot(real(gamma(1,:)),imag(gamma(1,:)),'b',...
    real(gamma(1,:)),-imag(gamma(1,:)),'b:',...
    real(gamma(2,:)),imag(gamma(2,:)),'r',...
    real(gamma(2,:)),-imag(gamma(2,:)),'r:',...
    'linewidth',2);
grid



