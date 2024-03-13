L1_nom=.001;
L1 = ureal('L1',L1_nom);
% L1 = L1_nom;
% ZL1=s*L1;
L2=.05;
% ZL2=s*L2;
L3=.01;
% ZL3=s*L3;
C1=.01;
% ZC1=1/(s*C1);
C2_nom=.005;
C2 = ureal('C2',C2_nom);
% C2 = C2_nom;
% ZC2=1/(s*C2);
C3=.001;
% ZC3=1/(s*C3);
R1=100;
% ZR1=R1;
R2=100;
% ZR2=R2;

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
num_states = size(G.A,1);
num_u_in = 2;
num_v_out = 2;
num_w_in = 4;
num_z_out = 3;
epsilon1 = 1;
epsilon2 = 1;

% Add Performance out and exogenous in
B1 = [zeros(num_states,num_w_in-1),G.B(:,1)];
C1 = [G.C(1,:);zeros(num_z_out-1,num_states)];
D11 = zeros(num_z_out,num_w_in);D11(1,1)=-1;
D12 = zeros(num_z_out,num_u_in);D12(2,1)=epsilon1;D12(3,2)=epsilon2;
D21 = [-1,1,0,0;0,0,1,0];
B = [B1,G.B];
C = [C1;G.C];
D = [D11,D12;D21,G.D];
P = ss(G.A,B,C,D);
% Add parametric uncertainty
[M,Delta] = lftdata(P);

% Parametric uncertainty without w or z
[M_b,Delta_b] = lftdata(G);
K = q_control(G_ss,eye(2));
cl = feedback(G*K,eye(2));
% Small Gain Theorem
norm(cl,inf)
norm(Delta_b.NominalValue,inf)

F = G.NominalValue*K;
Fss=ss(F);
syms k1 k2 real
s=sym('s');
n = length(Fss.a);
charpoly=simplify(det(s*eye(n)-(Fss.a-Fss.b*[k1 0;0 k2]*Fss.c)));
charpoly=collect(charpoly,s);
% pretty(charpoly);
charpoly=collect(charpoly,s);
k1_max = 1.2;
k1_min = .8;
k2_max = 1.2;
k2_min = .8;
charpoly1 = subs(charpoly,k1,k1_max);
charpoly1 = subs(charpoly1,k2,k2_max);
coeffs1 = coeffs(charpoly1,'All');
charpoly2 = subs(charpoly,k1,k1_min);
charpoly2 = subs(charpoly2,k2,k2_max);
coeffs2 = coeffs(charpoly2,'All');
charpoly3 = subs(charpoly,k1,k1_max);
charpoly3 = subs(charpoly3,k2,k2_min);
coeffs3 = coeffs(charpoly3,'All');
charpoly4 = subs(charpoly,k1,k1_max);
charpoly4 = subs(charpoly4,k2,k2_max);
coeffs4 = coeffs(charpoly4,'All');
ap = [];
an = [];

for i = 1:15
    ap(i) = max([coeffs1(i),coeffs2(i),coeffs3(i),coeffs4(i)]);
    an(i) = min([coeffs1(i),coeffs2(i),coeffs3(i),coeffs4(i)]);
end

poly1 = [];
poly2 = [];
poly3 = [];
poly4 = [];

for i = 1:15
    if mod(i,4)==0
        poly1(i) = an(i);
        poly2(i) = ap(i);
        poly3(i) = ap(i);
        poly4(i) = an(i);
    elseif mod(i,4)==1
        poly1(i) = an(i);
        poly2(i) = ap(i);
        poly3(i) = an(i);
        poly4(i) = ap(i);
    elseif mod(i,4)==2
        poly1(i) = ap(i);
        poly2(i) = an(i);
        poly3(i) = an(i);
        poly4(i) = ap(i);
    elseif mod(i,4)==3
        poly1(i) = ap(i);
        poly2(i) = an(i);
        poly3(i) = ap(i);
        poly4(i) = an(i);
    end
end
disp([max(real(roots(poly1))), max(real(roots(poly2))),max(real(roots(poly3))), max(real(roots(poly4)))]);

