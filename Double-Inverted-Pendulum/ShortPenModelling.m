%% Short Pendulum modelling
clc;
% clear;
syms m L l th dth ddth g Bp1 J Cp1 sgn1

% J = 1/12*m*L^2;

q = [th];
dq = [dth];
ddq = [ddth];

rp = [l*sin(th);-l*cos(th);0];
drp = jacobian(rp,q)*dq;

wp = [0;0;dth];

Ttot = simplify((1/2)*m*transpose(drp)*drp + (1/2)*transpose(wp)*J*wp);
Vtot = m*g*transpose(rp)*[0;1;0];

%% Define Mass Matrix
M = hessian(Ttot,dq);

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
G = simplify(G);

%% Define Coriolis Matrix
C = dM*dq - transpose(jacobian(Ttot,q));
C = simplify(C);

manipeqn = M*ddq + C + transpose(G) == -Bp1*dth

ddq_sol1 = simplify(M\(-Bp1*dth-C-transpose(G)))
%%
% L = 0.2096;
% l = 0.1143;
g = 9.81;

ddq_sol = subs(simplify(M\(-Bp1*dth-C-transpose(G))))