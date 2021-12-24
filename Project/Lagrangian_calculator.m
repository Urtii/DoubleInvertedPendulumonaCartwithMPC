syms x(t) theta(t) Q(t) m1 m2 m3 k l g ddx ddtheta ddQ I2 I3%Define the symbolic variables.
dx = diff(x,t);
%ddx = diff(dx,t);

dtheta = diff(theta,t);
%ddtheta = diff(dtheta,t);

dQ = diff(Q,t);
%ddQ = diff(dQ,t);

% g = 9.81;
% k = 0.5;
% l = 0.25;
% m1 = 0.6;
% m2 = 0.15;
% m3 = 0.3;

%Potential Energy
V = g * (k*cos(theta)*(m2+m3) - m3*l*cos(Q));
%Kinetic Energy
Kx = 1/2 * (m1*dx^2 + m2*(dx-k*dtheta*cos(theta))^2 + m3*(dx-k*dtheta*cos(theta)+l*cos(Q)*dQ));
Ky = 1/2 * (m2*(-k*dtheta*sin(theta))^2 + m3*(-k*dtheta*sin(theta)+l*sin(Q)*dQ)^2);
Kr = 1/2 * (I2*dtheta^2 + I3*dQ^2);

%Lagrangian Equation
L = simplify(expand(Kx + Ky + Kr - V));

%Take Lagrangian derivatives for all 3 variavles
fx = simplify(expand(diff(diff(L,dx),t) - diff(L,x)));
ftheta = simplify(expand(diff(diff(L,dtheta),t) - diff(L,theta)));
fQ = simplify(expand(diff(diff(L,dQ),t) - diff(L,Q)));

%Change diff(x(t),t,t) variables with ddx or corresponding variables
Fx = subs(subs(subs(fx,diff(x(t), t, t),ddx),diff(theta(t), t, t),ddtheta),diff(Q(t), t, t),ddQ);
Ft = subs(subs(subs(ftheta,diff(x(t), t, t),ddx),diff(theta(t), t, t),ddtheta),diff(Q(t), t, t),ddQ);
FQ = subs(subs(subs(fQ,diff(x(t), t, t),ddx),diff(theta(t), t, t),ddtheta),diff(Q(t), t, t),ddQ);

%Get solutions of ddx, ddtheta and ddQ in terms of position and velocity
eqns = [Fx, Ft, FQ];
vars = [ddx, ddtheta, ddQ];
S = solve(eqns,vars);
solx = simplify(S.ddx);
soltheta = simplify(S.ddtheta);
solQ = simplify(S.ddQ);

%Linearization
% linearize around (x=0,theta=0,Q=0)
Lin_ddx = subs(subs(subs(solx,x,0),theta,0),Q,0) + ...
    x*subs(subs(subs(diff(solx,x),x,0),theta,0),Q,0) + ...
    theta*subs(subs(subs(diff(solx,theta),x,0),theta,0),Q,0) + ...
    (Q)*subs(subs(subs(diff(solx,Q),x,0),theta,0),Q,0);

Lin_ddtheta = subs(subs(subs(soltheta,x,0),theta,0),Q,0) + ...
    x*subs(subs(subs(diff(soltheta,x),x,0),theta,0),Q,0) + ...
    theta*subs(subs(subs(diff(soltheta,theta),x,0),theta,0),Q,0) + ...
    (Q)*subs(subs(subs(diff(soltheta,Q),x,0),theta,0),Q,0);

Lin_ddQ = subs(subs(subs(solQ,x,0),theta,0),Q,0) + ...
    x*subs(subs(subs(diff(solQ,x),x,0),theta,0),Q,0) + ...
    theta*subs(subs(subs(diff(solQ,theta),x,0),theta,0),Q,0) + ...
    (Q)*subs(subs(subs(diff(solQ,Q),x,0),theta,0),Q,0);


A = [subs(Lin_ddx,)]