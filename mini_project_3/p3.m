%% Plant Declaration
L1=.001;
L2=.05;
L3=.01;
C1=.01;
C2=.005;
C3=.001;
R1=100;
R2=100;

A = [0,0,0,-1/L1,-1/L1,0;
     0,-R2/L2,0,0,1/L2,0;
     0,0,0,0,1/L3,-1/L3;
     1/C1,0,0,0,0,0;
     1/C2,-1/C2,-1/C2,0,-1/(C2*R1),0;
     0,0,1/C3,0,0,0];
B = [1/L1,0;
     0,-1/L2;
     0,0;
     0,0;
     1/(C2*R1),0;
     0,0];
C = [0,1,0,0,0,0;
     0,0,0,0,-1/R1,0];
D = [0,0;
     1/R1,0];

u_unc = [ureal('u1',0,'PlusMinus',.5),0;0,ureal('u2',0,'PlusMinus',.5)];

B = B*(eye(2)+u_unc);
D = D*(eye(2)+u_unc);

G = ss(A,B,C,D);



eps1=.0001;eps2=.0001;
n=length(G.a);
A=G.a;B2=G.b;C2=G.c;D22=G.d;
B1=[zeros(n,3) B2.NominalValue];
C1=[C2(1,:);zeros(2,n)];
D11=[-1 zeros(1,4);zeros(2,5)];
D12=[zeros(1,2);diag([eps1,eps2])];
D21=[-1 1 0 0 0;0 0 1 0 0];
P=ss(A,[B1 B2],[C1;C2],[D11 D12;D21 D22]);
nw=5;nu=2;nz=3;ny=2;

%% RGA (1)
f = freqresp(G,5);
rga = f.*pinv(f).'

%% 2a
% Wr=tf(5,[1 10])*tf(10,[1 .1]);
Wr=tf(5,[1 5])*tf(100,[1 .1]);
Wd=tf([2,0],[1,2,4]);
Wu=tf(5,[1 5])*.001;
% b1=20;b2=50;Wn=.1*[tf([1 0],[1 b1]) 0 ;0 tf([1 0],[1 b2])];
b1=15;b2=15;Wn=[tf([1 0],[1 b1]) 0 ;0 tf([1 0],[1 b2])];

% P(1,1) is the tracking error (from r to y1-r)
P_w = P;
P_w(2:3,:)=P(2:3,:)*Wu;
P_w(:,1)=P(:,1)*Wr;
P_w(:,2:3)=P(:,2:3)*Wn;
P_w(:,4)=P(:,4)*Wd;
P_w = minreal(P_w);
[K_inf,Twz_inf,gam_inf,infoinf]=hinfsyn(P_w,ny,nu);
[K_2,Twz_2,gam_2,info2]=h2syn(P_w,ny,nu);
Gcl_inf=lft(P,K_inf,nu,ny);
Gcl_2=lft(P,K_2,nu,ny);
figure(1);step(Gcl_inf(1,1),Gcl_2(1,1));legend('H_\infty','H_2');

n_Gcl_2_2 = norm(Gcl_2);
%% Simulink Plots (2a)
Kopt = -K_2;
n_max=[.1;.1];n_power=[.01;.01];
a1=10;a2=5;w1=5;w2=0;
Wn=Wn*eye(2);

figure(2);plot(out.tout,out.n,'linewidth',2);legend('n_1','n_2');
figure(3);plot(out.tout,out.e,'linewidth',2);legend('e_1','e_2');
figure(4);plot(out.tout,out.u,'linewidth',2);legend('u_1','u_2');
figure(5);plot(out.tout,out.up,'linewidth',2);legend('u_{p1}','u_{p2}');
figure(6);plot(out.tout,out.y,'linewidth',2);legend('y_1','y_2');


%% 2b
stable_inf = isstable(Gcl_inf);
stable_2 = isstable(Gcl_2);
[stab_marg_inf,wcu_inf] = robstab(Gcl_inf);
[stab_marg_2,wcu_2] = robstab(Gcl_2);

[Gcl_2b_inf,Delta_inf,blk_inf] = lftdata(Gcl_inf);
[Gcl_2b_2,Delta_2,blk_2] = lftdata(Gcl_2);
Gcl_2b_inf = minreal(Gcl_2b_inf);
Gcl_2b_2 = minreal(Gcl_2b_2);
figure(7);nyquist(Gcl_2b_inf(3,1));
figure(8);nyquist(Gcl_2b_inf(3,2));
figure(9);nyquist(Gcl_2b_2(3,1));
figure(10);nyquist(Gcl_2b_2(3,2));

bounds_inf = mussv(Gcl_2b_inf(1:2,1:2),[-1 0;-1 0]);
bounds_2 = mussv(Gcl_2b_2(1:2,1:2),[-1 0;-1 0]);

w = logspace(-5,5,100);

gamma=mvar_nyquist(tf(feedback(G*K_inf,eye(2))),w);
figure(11); 
plot(real(gamma(1,:)),imag(gamma(1,:)),'b',...
    real(gamma(1,:)),-imag(gamma(1,:)),'b:',...
    real(gamma(2,:)),imag(gamma(2,:)),'r',...
    real(gamma(2,:)),-imag(gamma(2,:)),'r:',...
    'linewidth',2);
grid

gamma=mvar_nyquist(feedback(G*K_2,eye(2)),w);
figure(12); 
plot(real(gamma(1,:)),imag(gamma(1,:)),'b',...
    real(gamma(1,:)),-imag(gamma(1,:)),'b:',...
    real(gamma(2,:)),imag(gamma(2,:)),'r',...
    real(gamma(2,:)),-imag(gamma(2,:)),'r:',...
    'linewidth',2);
grid
%% 3a
[M,Delta,BlockStructure] = lftdata(P);
% Wr=tf(5,[1 5])*tf(100,[1 .1]);
Wr = tf(.01,[1 .1]);
Wd=tf([2,0],[1,2,4]);
b1=15;b2=15;Wn=[tf([1 0],[1 b1]) 0 ;0 tf([1 0],[1 b2])];
M_w = M;
M_w(:,3) = M_w(:,3)*Wr;
M_w(:,4:5)=M(:,4:5)*Wn;
M_w(:,6)=M(:,6)*Wd;
M_w(:,7)=M(:,7)*Wd;

M_w = minreal(M_w);

[K_inf_unc,Twz_inf_unc,gam_inf_unc,infoinf_unc]=hinfsyn(M_w,ny,nu);
Gcl_3a=lft(P,K_inf_unc,nu,ny);
figure(13);step(Gcl_3a(1,1));
n_2 = norm(Gcl_3a(1,1));

Gcl_3a_unc = lftdata(Gcl_3a);
bounds_3a = mussv(Gcl_3a_unc(1:2,1:2),[-1 0;-1 0]);

%% 3b

P_w = P;
% Wr=20*tf(5,[1 5])*tf(.001,[1 .1]);
Wr = tf(.01,[1 .1]);
% Wr = tf(.01,[1 5]);
% Wd=tf([2,0],[1,2,4]);
% Wu=tf(5,[1 5])*.001;
% b1=15;b2=15;Wn=[tf([1 0],[1 b1]) 0 ;0 tf([1 0],[1 b2])];
% P_w(:,2:3)=P(:,2:3)*Wn;
% P_w(:,4)=P(:,4)*Wd;
P_w(:,1)=P(:,1)*Wr;
% P_w(2:3,:)=P(2:3,:)*Wu;
P_w = minreal(P_w);
[K_3b,CLperf] = musyn(P_w,ny,nu);
Gcl_3b=lft(P_w,K_3b);

figure(14);step(Gcl_3b(1,1));
Gcl_3b_unc = lftdata(Gcl_3b);
n_3b_inf_unc = norm(Gcl_3b_unc(1:2,1:2),'inf')
n_3b_inf = norm(Gcl_3b,'inf');
n_3b_2 = norm(Gcl_3b);

%% 4

Wr=tf(5,[1 5])*tf(100,[1 .01]);
Wd=tf([2,0],[1,2,4]);
Wu=tf(5,[1 5])*.001;
b1=15;b2=15;Wn=[tf([1 0],[1 b1]) 0 ;0 tf([1 0],[1 b2])];

% y1,u2
A = G(1,2).A;
B = G(1,2).B;
C = G(1,2).C;
D = G(1,2).D;
eps1=.0001;eps2=.0001;
n=length(G.a);
A=A;B2=B;C2=C;D22=D;
B1=zeros(n,5);
C1=[C2(1,:);zeros(2,n)];
D11=[-1 zeros(1,4);zeros(2,5)];
D12=[0;eps1;0];
D21=[-1 1 0 0 0];
P_y1_u2=ss(A,[B1 B2],[C1;C2],[D11 D12;D21 D22]);
P_w=P_y1_u2;
P_w(:,1)=P_y1_u2(:,1)*Wr;
P_w = minreal(P_w);

[K_y1_u2,Twz_2,gam_2a,info2]=h2syn(P_w,1,1);

Gcl1 = lft(P_y1_u2,K_y1_u2,1,1);
figure(15);
step(Gcl1(1,1));

% y2,u1
A = G(2,1).A;
B = G(2,1).B;
C = G(2,1).C;
D = G(2,1).D;
eps1=.0001;eps2=.0001;
n=length(G.a);
A=A;B2=B;C2=C;D22=D;
B1=zeros(n,5);
C1=[zeros(3,n)];
D11=[-1 zeros(1,4);zeros(2,5)];
D12=[0;0;eps2];
D21=[0 0 1 0 0];
P_y2_u1=ss(A,[B1 B2],[C1;C2],[D11 D12;D21 D22]);
P_w=P_y2_u1;
P_w(:,1)=P_y2_u1(:,1)*Wr;
P_w = minreal(P_w);

[K_y2_u1,Twz_2,gam_2b,info2]=h2syn(P_w,1,1);
Gcl2 = lft(P_y2_u1,K_y2_u1,1,1);

cl1=feedback(P,K_y1_u2,7,4,+1);
cl_siso = feedback(cl1,K_y2_u1,6,5,+1);
figure(16);
step(cl_siso(1,1));

