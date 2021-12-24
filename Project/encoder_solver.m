syms x(t) theta(t) Q(t) m1 m2 m3 k l g ddx ddtheta ddQ

dx = diff(x,t);
dtheta = diff(theta,t);
dQ = diff(Q,t);

Fx = 0 == m1*ddx + m2*ddx + k*m2*sin(theta(t))*diff(theta(t), t)^2 - k*m2*cos(theta(t))*ddtheta;
Ft = 0 == k^2*m2*ddtheta + k^2*m3*ddtheta + (k^2*m3*sin(2*theta(t))*diff(theta(t), t)^2)/2 - k^2*m3*cos(theta(t))^2*ddtheta - k*m2*cos(theta(t))*ddx - g*k*m2*sin(theta(t)) - g*k*m3*sin(theta(t)) - k*l*m3*cos(Q(t))*sin(theta(t))*diff(Q(t), t)^2 - k*l*m3*sin(Q(t))*sin(theta(t))*ddQ;
FQ = 0 == - m3*l^2*cos(Q(t))^2*ddQ + m3*l^2*ddQ - k*m3*sin(Q(t))*sin(theta(t))*l*ddtheta + m3*sin(Q(t))*l^2*cos(Q(t))*diff(Q(t), t)^2 - k*m3*sin(Q(t))*cos(theta(t))*l*diff(theta(t), t)^2 + g*m3*sin(Q(t))*l;

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
    Q*subs(subs(subs(diff(solx,Q),x,0),theta,0),Q,0);

Lin_ddtheta = subs(subs(subs(soltheta,x,0),theta,0),Q,0) + ...
    x*subs(subs(subs(diff(soltheta,x),x,0),theta,0),Q,0) + ...
    theta*subs(subs(subs(diff(soltheta,theta),x,0),theta,0),Q,0) + ...
    Q*subs(subs(subs(diff(soltheta,Q),x,0),theta,0),Q,0);

% Lin_ddQ = subs(subs(subs(solQ,x,0),theta,0),Q,0) + ...
%     x*subs(subs(subs(diff(solQ,x),x,0),theta,0),Q,0) + ...
%     theta*subs(subs(subs(diff(solQ,theta),x,0),theta,0),Q,0) + ...
%     Q*subs(subs(subs(diff(solQ,Q),x,0),theta,0),Q,0);