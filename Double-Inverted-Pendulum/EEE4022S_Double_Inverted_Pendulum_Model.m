%% Benjamin Measures
% EEE4022S Project
% Linear Double Inverted Pendulum

%% Modeling the system
% Co-ordinates
clc;
clear;
syms x dx ddx th1 dth1 ddth1 th2 dth2 ddth2
q = [th1;th2;x];
dq = [dth1;dth2;dx];
ddq = [ddth1;ddth2;ddx];

% Parameters
syms mc m1 m2 l1 l2 L1 L2

syms Fc Bc Bp1 Bp2 Beq J1 J2
B_damp = [Bp1,0,0;0,Bp2,0;0,0,Bc];

syms R_m k_t n_m k_m K_g n_g r_mp V_m

syms cp1 cp2 sgn1 sgn2

% Constants
syms g

% Rotations
R01 = Rotz(th1);
R10 = transpose(R01);
R12 = Rotz(th2);
R02 = R12*R01;
R20 = transpose(R02);

%% Mass Positions
% Cart
rC_0 = [x;0;0];

% Pendulum 1
rP1_1 = [0;-l1;0];
rP1_0 = simplify(rC_0 + R10*rP1_1);

% Pendulum 2
rP2_2 = [0;-l2;0];
rP2_0 = simplify(rC_0 + R10*[0;-L1;0] + R20*rP2_2);

%% Mass Velocities
drC = jacobian(rC_0,q)*dq;
drP1 = jacobian(rP1_0,q)*dq;
drP2 = simplify(expand(jacobian(rP2_0,q)*dq));

%% Angular Velocity
w01_1 = [0;0;dth1];
w01_2 = R12*w01_1;

w12_2 = [0;0;dth2];
w02_2 = w01_2 + w12_2;

%% Kinetic Energy
TC = (1/2)*mc*transpose(drC)*drC;
TP1 = (1/2)*m1*transpose(drP1)*drP1 + (1/2)*transpose(w01_1)*J1*w01_1;
TP2 = (1/2)*m2*transpose(drP2)*drP2 + (1/2)*transpose(w02_2)*J2*w02_2;
Ttot = simplify(TC+TP1+TP2);

%% Potential Energy
VP1 = m1*g*transpose(rP1_0)*[0;1;0];
VP2 = m2*g*transpose(rP2_0)*[0;1;0];
Vtot = simplify(VP1+VP2);

%% Define Mass Matrix
M = simplify(hessian(Ttot,dq));
%% Define Mass Matrix Deriv
dM = sym(zeros(length(M),length(M)));
for i=1:length(M)
    for j=1:length(M)
        dM(i,j) = jacobian(M(i,j),q)*dq;
    end
end
dM = simplify(dM);

%% Define Gravity Matrix
G = jacobian(Vtot,q);
G = simplify(transpose(G));

%% Define Coriolis Matrix
C = dM*dq - transpose(jacobian(Ttot,q));
C = simplify(C);

%%
M11(1,1) = M(1,1);
M11(1,2) = M(1,2);
M11(2,1) = M(2,1);
M11(2,2) = M(2,2);

N1(1,1) = C(1,1) + G(1,1);
N1(2,1) = C(2,1) + G(2,1);

N2 = C(3,1) + G(3,1);

M12(1,1) = M(1,3);
M12(2,1) = M(2,3);

M21 = transpose(M12);

M22 = M(3,3);

M22_bar = simplify(M22 - M21*(M11\M12));
N2_bar = N2-M21*(M11\N1);

ddq1 = [ddth1;ddth2];
ddq2 = ddx;

Ep = simplify(subs(TP1+TP2+Vtot,[x dx],[0 0]));
Ep1 = simplify((1/2)*transpose(w01_1)*J1*w01_1 + (1/2)*transpose(w02_2)*J2*w02_2 + m1*g*transpose(rP1_0)*[0;1;0] + m2*g*transpose(rP2_0)*[0;1;0]);

Er = subs(TP1+TP2+Vtot,[th1 th2 x dth1 dth2 dx],[0 0 0 0 0 0]);
syms u f1

g_ksi =  simplify([0; 0; -M11\M12]);
u_bar = simplify(-[jacobian(Ep,th1),jacobian(Ep,th2),jacobian(Ep,dth1),jacobian(Ep,dth2)]*g_ksi);

V_Fc_relationship = f1==(n_g*K_g*n_m*k_t*(V_m*r_mp-K_g*k_m*dx))/(R_m*(r_mp^2));
V_sol = expand(solve(V_Fc_relationship,V_m));
V_sol = subs(V_sol,f1,mc*u);
%%
E_test = simplify(subs(TP1 + TP2 + Vtot,[x dx],[0 0]));
E_test1 = subs(E_test,[dth1 dth2],[0 0]);
E_test2 = simplify(subs(E_test,[th1 th2],[pi/2 0]));
%% Manipulator Equation
ManipulatorEqn = simplify(M*ddq + C + transpose(G)) == [0;0;Fc] - B_damp*dq;
ddq_final = simplify(M\([0;0;Fc] - B_damp*dq - C - transpose(G)));

%% Sub values in
Fc = (n_g*K_g*n_m*k_t*(V_m*r_mp-K_g*k_m*dx))/(R_m*(r_mp^2));
g = 9.81;
l1 = 0.1483;
L1 = 0.2096;
l2 = 0.1778;
m1 = 0.25244;
m2 = 0.10562;
mc = 0.57;
J1 = 0.001;
J2 = 8.2572e-04;
Bp1 = 5.6134e-04;
Bp2 = 4.2362e-04;
Cp1 = 0;
Cp2 = 0;
sgn1 = 0;
sgn2 = 0;
Bc = 9.93378;

R_m = 2.6; 
k_t = 0.00767;
n_m = 1;
k_m = 0.00767;
K_g = 3.71;
n_g = 1;
r_mp = 6.35*10^(-3);

ddq_final = M\([0;0;Fc]-B_damp*dq-C-G);
ddq_final = simplify(subs(ddq_final));

%% State Space
dyn = [dth1;dth2;dx;ddq_final(1);ddq_final(2);ddq_final(3)];
A = jacobian(dyn,[th1;th2;x;dth1;dth2;dx]);
A = subs(A,[th1,th2,x,dth1,dth2,dx,V_m],[pi,0,0,0,0,0,0]);
A = double(A);

B = jacobian(dyn,V_m);
B = subs(B,[th1,th2,x,dth1,dth2,dx,V_m],[pi,0,0,0,0,0,0]);
B = double(B);

C1 = [
    1 0 0 0 0 0
    0 1 0 0 0 0
    0 0 1 0 0 0
    ];

D = [
    0
    0
    0
    ];

sys = ss(A,B,C1,D);

%%
Q = eye(6);
Q(1,1) = 60;
Q(2,2) = 60;
Q(3,3) = 10;
Q(4,4) = 0.1;
Q(5,5) = 0.1;
Q(6,6) = 0;

R = 0.01;

K = lqr(A,B,Q,R)

wcf_1 = 2 * pi * 50.0;    
zetaf_1 = 0.9;          
wcf_2 = 2 * pi * 10;     
zetaf_2 = 0.9;          

%%
T = [B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B];
rank(T)

obsv = [C1;C1*A;C1*(A^2);C1*(A^3);C1*(A^4);C1*(A^5)];
rank(obsv)  
%% Rotation Matrices
function A = Rotz(th)
    A = [cos(th)   sin(th) 0;... 
         -sin(th)  cos(th) 0;...
         0        0        1];
end