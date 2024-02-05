% ========================================================================
% PLANT:
% ========================================================================
p = @(varargin)1/sum(1./[varargin{:}]);
syms s;

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

Z1=p(ZC1+ZL1,ZR1);
Z2=ZR2+ZL2;
Z3=p(ZL3+ZC3,ZC2);
Ztot = p(Z2,Z3)/Z2;
G = Ztot/(p(Z2,Z3)+Z1);
[Gnumsym,Gdensym] = numden(G);
Gnum=sym2poly(Gnumsym);
Gden=sym2poly(Gdensym);
plant=tf(Gnum,Gden)
[A,B,C,D] = tf2ss(Gnum,Gden);
% figure()
% bode(plant)
% figure()
% step(plant)
plant_data = stepinfo(plant);
% figure()
% rlocus(plant)

% ========================================================================
% STEP RESPONSE CONTROLLER:
% ========================================================================
resonant_peak = 566;
a = 30;
notch = 25*tf([1,0,resonant_peak^2],resonant_peak^2)*tf([a*resonant_peak],[1,a*resonant_peak])*tf([resonant_peak/a],[1,resonant_peak/a])
pid1 = pid(500,300,0);
ol_controller1 = notch*pid1
open_loop1 = ol_controller1*plant;
figure(1)
bode(notch)
figure(2)
bode(pid1)
figure(3)
margin(open_loop1)
output1 = feedback(open_loop1,1)
% pidtune(notch*plant,'PID')
figure(4)
step(output1)
controlled_data=stepinfo(output1);
figure(5)
bode(output1)

[n,d]=tfdata(output1,'v');
[step_controller_A,step_controller_B,step_controller_C,step_controller_D] = tf2ss(n,d);
[step_controller_z,step_controller_p,step_controller_k] = tf2zp(n,d);
Gdy=feedback(plant,ol_controller1);
Gne=feedback(-1,plant*ol_controller1*-1);
figure(6)
bode(Gdy)
figure(7)
bode(Gne)

% ========================================================================


% ========================================================================
% TRACKING CONTROLLER:
% ========================================================================

noise_peak = 20;
a = 1;
notch_noise = .5*tf([1,0,noise_peak^2],noise_peak^2)*tf([a*noise_peak],[1,a*noise_peak])*tf([noise_peak/a],[1,noise_peak/a]);
ol_controller2 = ol_controller1*notch_noise;
figure(8)
open_loop2 = ol_controller2*plant;
margin(open_loop2)
figure(9)
bode(feedback(open_loop2,1))


