%% Load the plant
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
G = ss(A,B,C,D);
eps1=.0001;eps2=.0001;
n=length(G.a);
A=G.a;B2=G.b;C2=G.c;D22=G.d;
B1=zeros(n,1);
C1=[C2(1,:);zeros(2,n)];
D11=[-1;0;0];
D12=[zeros(1,2);diag([eps1,eps2])];
D21=[-1;0];
P_truth=ss(A,[B1 B2],[C1;C2],[D11 D12;D21 D22]);

%% Sytem ID
ts=.001;T=50;Ns=fix(T/ts); % 1ms sampling period
fmax=1/ts*(Ns/2-1)/Ns;fs=1/ts/Ns; % 1/2 of the frequency response

fband=500;
Nband=fix(fband/fs);

uf_mag=zeros(1,Ns/2);
uf_mag(1:Nband)=1;

usat=10;
[u_schroed,t,mags,phs]=schroed(ts,Ns,uf_mag,usat);
tvec=(0:Ns-1)*ts;
figure(1),plot(tvec,u_schroed);
title('Schroder phase signal (time domain)');
xlabel('time (sec)');ylabel('input (N)')

fvec=(0:Ns/2-1)*fs;
fvec_band=fvec(1:Nband);
u_f=t2f(u_schroed)*ts;
figure(2),plot(fvec,abs(u_f),'x');
axis([0 2*fband 0 5])
title('Schroder phase signal (frequency domain)');
xlabel('frequency (Hz)');ylabel('input magnitude');
figure(3),plot(fvec,unwrap(angle(u_f)),'o');
axis([0 2*fband 0 200])
title('Schroder phase signal (frequency domain)');
xlabel('frequency (Hz)');ylabel('input phase');

numrun=10;
u_std=1;
y_std=0.1;
u_excite=kron(ones(numrun,1),u_schroed')';

tvec_numrun=(0:numrun*Ns-1)*ts;
u_noise1=randn(1,length(tvec)*numrun)*u_std;
y_noise1=randn(1,length(tvec)*numrun)*y_std;
y_u1_full=lsim(G(:,1),u_excite+u_noise1,tvec_numrun)+y_noise1';
u_noise2=randn(1,length(tvec)*numrun)*u_std;
y_noise2=randn(1,length(tvec)*numrun)*y_std;
y_u2_full=lsim(G(:,2),u_excite+u_noise2,tvec_numrun)+y_noise2';

y1_u1=mean(reshape(y_u1_full(:,1),Ns,numrun)');
y2_u1=mean(reshape(y_u1_full(:,2),Ns,numrun)');
y1_u2=mean(reshape(y_u2_full(:,1),Ns,numrun)');
y2_u2=mean(reshape(y_u2_full(:,2),Ns,numrun)');

figure(4),plot(tvec,y1_u1,tvec,y1_u2,tvec,y2_u1,tvec,y2_u2);

y11=t2f(y1_u1)*ts/abs(u_f(2));
y21=t2f(y2_u1)*ts/abs(u_f(2));
y12=t2f(y1_u2)*ts/abs(u_f(2));
y22=t2f(y2_u2)*ts/abs(u_f(2));

[mag,phase]=bode(G,2*pi*fvec);

figure(5),semilogx(fvec,squeeze(mag(1,1,:)),'o',fvec,abs(y11),'x');
figure(6),semilogx(fvec,squeeze(mag(2,1,:)),'o',fvec,abs(y21),'x');
figure(7),semilogx(fvec,squeeze(mag(1,2,:)),'o',fvec,abs(y12),'x');
figure(8),semilogx(fvec,squeeze(mag(2,2,:)),'o',fvec,abs(y22),'x');

%% Freq Domain ID

fid_band=250;
Nid_band=fix(fid_band/fs);
G_freq=zeros(2,2,Nid_band);
G_freq(1,1,:)=y11(1:Nid_band);
G_freq(1,2,:)=y12(1:Nid_band);
G_freq(2,1,:)=y21(1:Nid_band);
G_freq(2,2,:)=y22(1:Nid_band);
G_frd=frd(G_freq,fvec(1:Nid_band),'frequencyunit','Hz');
figure(9),bode(G,G_frd)
dat_f=iddata(G_frd,'Domain','Frequency');
% G_est_f=ssest(dat_f,6);
G_est=n4sid(dat_f,10);
% figure(10);compare(dat_f,G_est_f);
figure(10);bode(G,G_est);legend('true','identified');
% present(G_est_f)
% % G_disc=ss(G_est_f.a,G_est_f.b,G_est_f.c,G_est_f.d,ts);
% G_cont_f=d2c(G_est_f,'tustin');
%% Time Domain ID
u_excite=kron(ones(numrun,1),u_schroed')';
u_in=[u_excite;u_excite];
y_noise1=randn(1,length(tvec)*numrun)*y_std;
y_noise2=randn(1,length(tvec)*numrun)*y_std;
y_full=lsim(G,u_in,tvec_numrun)+[y_noise1' y_noise2'];
y1=mean(reshape(y_full(:,1),Ns,numrun)');
y2=mean(reshape(y_full(:,2),Ns,numrun)');
dat=iddata([y1' y2'],[u_schroed' u_schroed'],ts);
G_est=ssest(dat,6,'Ts',ts);
present(G_est)
figure(11),compare(dat,G_est);
figure(12);bode(G,G_est);legend('true','id');

%% Determine Range of Uncertainty in model
u_std=1;
y_std=0.1;
inf_norms = [];
u_excite=kron(ones(numrun,1),u_schroed')';
u_in=[u_excite;u_excite];
tvec_numrun=(0:numrun*Ns-1)*ts;
for i  =1:50
    y_noise1=randn(1,length(tvec)*numrun)*y_std;
    y_noise2=randn(1,length(tvec)*numrun)*y_std;
    y_full=lsim(G,u_in,tvec_numrun)+[y_noise1' y_noise2'];
    y1=mean(reshape(y_full(:,1),Ns,numrun)');
    y2=mean(reshape(y_full(:,2),Ns,numrun)');
    dat=iddata([y1' y2'],[u_schroed' u_schroed'],ts);
    % G_est=ssest(dat,6,'Ts',ts);
    G_est=ssest(dat,6,'Ts',ts);
    present(G_est)
%     figure(13),compare(dat,G_est);
%     G_disc=ss(G_est.a,G_est.b,G_est.c,G_est.d,ts);
%     G_cont=d2c(G_disc,'tustin');
    inf_norms(i)=norm(G_est,'inf');
end
mean_gain = mean(inf_norms);
std_gain = std(inf_norms);
range_gain = max(inf_norms)-min(inf_norms);
%% Create Uncertain Object
% Load best model
load working.mat
out_unc = [ureal('y1',0)*range_gain,0;0,ureal('y2',0)*range_gain];
G_combined = G_est*(eye(2)+out_unc);
bode(G_combined);

eps1=.001;eps2=.001;
n=length(G_combined.a);
A=G_combined.a;B2=G_combined.b;C2=G_combined.c;D22=G_combined.d;
B1=zeros(n,1);
C1=[C2(1,:);zeros(2,n)];
D11=[-1;0;0];
D12=[zeros(1,2);diag([eps1,eps2])];
D21=[-1;0];
P=ss(A,[B1 B2],[C1;C2],[D11 D12;D21 D22]);

%% Create Controller
P_w = P;
Wr = tf(5,[1 5])*tf(1,[1 .01]);
P_w(1,:)=P(1,:)*Wr;
P_w = minreal(P_w);
[K,CLperf] = musyn(P_w,2,2);

Gcl = lft(P_truth,K);
figure(13);
step(Gcl(1,1));

%% Reduce order

[Kb,S]=balreal(K);
elim = (S<200);
Kr = xelim(Kb,6:8,'truncate');

Gcl_reduce = lft(P_truth,Kr);
figure(14);
step(Gcl_reduce(1,1));
Gcl_r_t = Gcl_reduce;
Gcl_r_t.D = 0
Gcl_r_t(:,1) = Gcl_r_t(:,1)*tf(10,[1 65]);
Gcl_t = Gcl;
Gcl_t.D = 0
Gcl_t(:,1) = Gcl_t(:,1)*tf(10,[1 65]);
ol_r_norm = norm(Gcl_r_t)
ol_norm = norm(Gcl_t)


