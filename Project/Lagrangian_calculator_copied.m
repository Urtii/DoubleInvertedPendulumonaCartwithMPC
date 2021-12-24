syms theta0(t) theta1(t) theta2(t) m0 m1 m2 L1 L2 I1 I2 g ...
    %Define the symbolic variables.
dtheta0 = diff(theta0,t);
dtheta1 = diff(theta1,t);
dtheta2 = diff(theta2,t);

g = 9.81;
L1 = 0.5;
L2 = 0.25;
m0 = 0.6;
m1 = 0.15;
m2 = 0.3;
I1 = 1/12*m1*L1^2;
I2 = 1/12*m2*L2^2;

%Kinetic Energy
K0 = 1/2*m0*dtheta0^2;

K1 = 1/2*m1*dtheta0^2 + 1/2*(m1*L1^2+I1)*dtheta1^2 + m1*L1*dtheta0*dtheta1*cos(theta1);

K2 = 1/2*m2*dtheta0^2 + 1/2*m2*L1^2*dtheta1^2 + 1/2*(m2*L2^2+I2)*dtheta2^2 + ...
    m2*L1*dtheta0*dtheta1*cos(theta1) + m2*L2*dtheta0*dtheta2*cos(theta2) + ...
    m2*L1*L2*dtheta1*dtheta2*cos(theta1+theta2);
%Potential Energy
V1 = g*m1*L1*cos(theta1);
V2 = g*m2*(L1*cos(theta1) - L2*cos(theta2));


%Lagrangian Equation
L = simplify(expand(K0+K1+K2 - V1-V2));

%Take Lagrangian derivatives for all 3 variavles
ftheta0 = simplify(expand(diff(diff(L,dtheta0),t) - diff(L,theta0)));
ftheta1 = simplify(expand(diff(diff(L,dtheta1),t) - diff(L,theta1)));
ftheta2 = simplify(expand(diff(diff(L,dtheta2),t) - diff(L,theta2)));


syms ddtheta0 ddtheta1 ddtheta2
%Change diff(theta0(t),t,t) variables with ddtheta0 or corresponding variables
Ftheta0 = subs(subs(subs(ftheta0,diff(theta0(t), t, t),ddtheta0),diff(theta1(t), t, t),ddtheta1),diff(theta2(t), t, t),ddtheta2);
Ftheta1 = subs(subs(subs(ftheta1,diff(theta0(t), t, t),ddtheta0),diff(theta1(t), t, t),ddtheta1),diff(theta2(t), t, t),ddtheta2);
Ftheta2 = subs(subs(subs(ftheta2,diff(theta0(t), t, t),ddtheta0),diff(theta1(t), t, t),ddtheta1),diff(theta2(t), t, t),ddtheta2);

Ftheta0 = subs(subs(subs(Ftheta0,diff(theta0(t), t),dtheta0),diff(theta1(t), t),dtheta1),diff(theta2(t), t),dtheta2);
Ftheta1 = subs(subs(subs(Ftheta1,diff(theta0(t), t),dtheta0),diff(theta1(t), t),dtheta1),diff(theta2(t), t),dtheta2);
Ftheta2 = subs(subs(subs(Ftheta2,diff(theta0(t), t),dtheta0),diff(theta1(t), t),dtheta1),diff(theta2(t), t),dtheta2);

%Get solutions of ddtheta0, ddtheta1 and ddtheta2 in terms of position and velocity
eqns = [Ftheta0, Ftheta1, Ftheta2];
vars = [ddtheta0, ddtheta1, ddtheta2];
S = solve(eqns,vars);
soltheta0 = simplify(S.ddtheta0);
soltheta1 = simplify(S.ddtheta1);
soltheta2 = simplify(S.ddtheta2);

%Linearization
% linearize around (theta0=0,theta1=0,theta2=0)
Ldd0 = simplify(subs(subs(subs(soltheta0,theta0,0),theta1,0),theta2,0) + ...
    theta0*subs(subs(subs(diff(soltheta0,theta0),theta0,0),theta1,0),theta2,0) + ...
    theta1*subs(subs(subs(diff(soltheta0,theta1),theta0,0),theta1,0),theta2,0) + ...
    (theta2)*subs(subs(subs(diff(soltheta0,theta2),theta0,0),theta1,0),theta2,0));

Ldd1 = simplify(subs(subs(subs(soltheta1,theta0,0),theta1,0),theta2,0) + ...
    theta0*subs(subs(subs(diff(soltheta1,theta0),theta0,0),theta1,0),theta2,0) + ...
    theta1*subs(subs(subs(diff(soltheta1,theta1),theta0,0),theta1,0),theta2,0) + ...
    (theta2)*subs(subs(subs(diff(soltheta1,theta2),theta0,0),theta1,0),theta2,0));

Ldd2 = simplify(subs(subs(subs(soltheta2,theta0,0),theta1,0),theta2,0) + ...
    theta0*subs(subs(subs(diff(soltheta2,theta0),theta0,0),theta1,0),theta2,0) + ...
    theta1*subs(subs(subs(diff(soltheta2,theta1),theta0,0),theta1,0),theta2,0) + ...
    (theta2)*subs(subs(subs(diff(soltheta2,theta2),theta0,0),theta1,0),theta2,0));

A= [0  1  0                                  0  0                                  0;
    0  0  subs(subs(Ldd0,theta1,1),theta2,0) 0  subs(subs(Ldd0,theta1,0),theta2,1) 0;
    0  0  0                                  1  0                                  0;
    0  0  subs(subs(Ldd1,theta1,1),theta2,0) 0  subs(subs(Ldd1,theta1,0),theta2,1) 0;
    0  0  0                                  0  0                                  1;
    0  0  subs(subs(Ldd2,theta1,1),theta2,0) 0  subs(subs(Ldd2,theta1,0),theta2,1) 0];

% B =[0 0;
%     5 7;
%     0 0;
%     3 11;
%     0 0;
%     2 13];
B = [0;3;0;5;0;7]
[B A*B A^2*B A^3*B A^4*B A^5*B]
det([B A*B A^2*B A^3*B A^4*B A^5*B])