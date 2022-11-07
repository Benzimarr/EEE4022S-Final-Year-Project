%% Benjamin Measures
% EEE4022S Project
% Linear Inverted Pendulum

%% Parameters
clc;
clear;

syms x th dx dth ddx ddth
q = [th;x];
dq = [dth;dx];
ddq = [ddth;ddx];

syms mc mp J Bc Bp g l L Fc

syms R_m k_t n_m k_m K_g n_g r_mp V_m

R01 = Rotz(th);
R10 = transpose(R01);

%% Mass positions
rC_0 = [x;0;0];

rP_1 = [0;-l;0];
rP_0 = rC_0 + R10*rP_1;

%% Velocities
drC = jacobian(rC_0,q)*dq;
drP = jacobian(rP_0,q)*dq;

%% Angular Velocity
w01_1 = [0;0;dth];

%% Kinetic Energy
TC = (1/2)*mc*transpose(drC)*drC;
TP = simplify((1/2)*mp*transpose(drP)*drP + (1/2)*transpose(w01_1)*J*w01_1);
Ttot = (TC+TP);

%% Potential Energy
Vtot = mp*g*[0,1,0]*rP_0;

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

%% Manipulator Dynamics
Fc = (n_g*K_g*n_m*k_t*(V_m*r_mp-K_g*k_m*dx))/(R_m*(r_mp^2));

Manip_Eqn = M*ddq + C + G == [-Bp*dth;Fc-Bc*dx]; 
ddq_sol = simplify(subs(M\([-Bp*dth;Fc-Bc*dx]-C-G)));

%% Energy Control
% Energy for pendulum
E1 = (1/2)*transpose(w01_1)*J*w01_1+Vtot;
E = subs(Ttot+Vtot,[x dx],[0 0]);
Er = subs(E,[th dth],[pi 0]);
E_dot = simplify(jacobian(E,q)*dq + jacobian(E,dq)*ddq);

eqn5 = simplify(subs(E_dot,ddth,eqn4));

%%
syms u Fc1
VFc_relationship = Fc1 ==(n_g*K_g*n_m*k_t*(V_m*r_mp-K_g*k_m*dx))/(R_m*(r_mp^2));
V_sol = expand(solve(VFc_relationship,V_m));
V_sol = subs(V_sol,Fc1,mc*u);

%%
g = 9.81;
mp = 0.230;
mc = 0.57;
l = 0.3302;
Bc = 9.3378;
Bp = 0.0019;
J = 0.0073;

R_m = 2.6; 
k_t = 0.00767;
n_m = 1;
k_m = 0.00767;
K_g = 3.71;
n_g = 1;
r_mp = 6.35*10^(-3);
%%
ddq_sol = simplify(subs(ddq_sol));

%%
dyn = [dth;dx;ddq_sol(1);ddq_sol(2)];
A = jacobian(dyn,[th;x;dth;dx]);
A = subs(A,[th,x,dth,dx,V_m],[pi,0,0,0,0]);
A = double(A);

B = jacobian(dyn,V_m);
B = subs(B,[th,x,dth,dx,V_m],[pi,0,0,0,0]);
B = double(B);

C1 = [
    1 0 0 0
    0 1 0 0
    ];

D = [
    0
    0
    ];

sys = ss(A,B,C1,D);

%%
Q = eye(4);
Q(1,1) = 60;
Q(2,2) = 20;
Q(3,3) = 0.1;
Q(4,4) = 0;

R = 0.01;

K = lqr(A,B,Q,R)

%%
T = [B A*B (A^2)*B (A^3)*B];
rank(T)

wcf = 2 * pi * 10.0;  
zetaf = 0.9;        
%% Rotation Matrices
function A = Rotz(th)
    A = [cos(th)   sin(th) 0;... 
         -sin(th)  cos(th) 0;...
         0        0        1];
end