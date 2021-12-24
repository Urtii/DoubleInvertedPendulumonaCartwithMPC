clc, clear all, close all;
%% Vehicle Model

fi = 0.03; 
tao = 0.05;
Ts = 0.01;
tend = 250;
N = 8;
% x = [x; theta; alpha; dx; dtheta; dalpha]
LAc =   [0,            0,            0, 1, 0, 0;
         0,            0,            0, 0, 1, 0;
         0,            0,            0, 0, 0, 1;
         0,  -26487/4900,   2943/12250, 0, 0, 0;
         0,   70632/1225,   35316/1225, 0, 0, 0;
         0, -105948/1225, -553284/6125, 0, 0, 0];

LBc =   [      0,        0;
               0,        0;
               0,        0;
         193/196,   97/392;
          -90/49, -1485/49;
          -12/49,  4506/49];

Cc = [1 0.5 0.5 0.001 0.1 0.1];
Dc = 0;
sys = ss(LAc,LBc,Cc,Dc);
sys_disc = c2d(sys,Ts);
[A,B,C,D] = ssdata(sys_disc);

x0 = [1;0;0;0;0;0];
% Optimal control solution for N = 8
G = [zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2);
     B          zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2);
     A*B        B          zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2);
     A^2*B      A*B        B          zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2);
     A^3*B      A^2*B      A*B        B          zeros(6,2) zeros(6,2) zeros(6,2) zeros(6,2);
     A^4*B      A^3*B      A^2*B      A*B        B          zeros(6,2) zeros(6,2) zeros(6,2);
     A^5*B      A^4*B      A^3*B      A^2*B      A*B        B          zeros(6,2) zeros(6,2);
     A^6*B      A^5*B      A^4*B      A^3*B      A^2*B      A*B        B          zeros(6,2);
     A^7*B      A^6*B      A^5*B      A^4*B      A^3*B      A^2*B      A*B        B];
H = [eye(6); A; A^2; A^3; A^4; A^5; A^6; A^7; A^8];
Q = C'*C;
R = [0.01 0;
     0    1];
%Q = eye(3);
Pinf = idare(A,B,Q,R,zeros(6,2),eye(6) );
Kinf = inv(R+B'*Pinf*B)*B'*Pinf*A;
% A*X*A' - X + Q = 0;  X = dlyap(A,Q)
P = dlyap( (A-B*Kinf)',Q+Kinf'*R*Kinf);
Qf = P;
Qbar = blkdiag(Q,Q,Q,Q,Q,Q,Q,Q,Qf);
Rbar = blkdiag(R,R,R,R,R,R,R,R);
M = G'*Qbar*G + Rbar;
% input bound: umin <= u <= umax
u1min = -3;
u1max = 3;
u2min = -4;
u2max = 4;
lb = [u1min;u2min;u1min;u2min;u1min;u2min;u1min;u2min;u1min;u2min;u1min;u2min;u1min;u2min;u1min;u2min];
ub = [u1max;u2max;u1max;u2max;u1max;u2max;u1max;u2max;u1max;u2max;u1max;u2max;u1max;u2max;u1max;u2max];
% Apply MPC steps
xVec(:,1) = x0;
yVec(1) = C*x0;
uVec = [0;0];
for kk = 1:250
    alpha = G'*Qbar'*H*xVec(:,kk);
    Usol = quadprog(M,alpha',[],[],[],[],lb,ub);
    uVec(:,kk) = [Usol(1);Usol(2)];
    xVec(:,kk+1) = A*xVec(:,kk) + B*uVec(:,kk);
    yVec(kk+1) = C*xVec(:,kk+1);
    Xsol(:,1) = xVec(:,kk);
    Xsol(:,2) = A*Xsol(:,1) + B*[Usol(1);Usol(2)];
    Xsol(:,3) = A*Xsol(:,2) + B*[Usol(1);Usol(2)];
    Xsol(:,4) = A*Xsol(:,3) + B*[Usol(1);Usol(2)];
    Ysol(1) = C*Xsol(:,1);
    Ysol(2) = C*Xsol(:,2);
    Ysol(3) = C*Xsol(:,3);
    Ysol(4) = C*Xsol(:,4);
end


uVec = [uVec uVec(:,end)];
tVec = [0:Ts:250*Ts];
% figure;
figure, subplot(3,1,1)
stairs(tVec,uVec(1,:),'LineWidth',2);
hold all;
xlabel('time [sec]')
grid
ylabel('u0')
title('Input u0')
subplot(3,1,2)
stairs(tVec,uVec(2,:),'LineWidth',2)
hold all;
grid
xlabel('time [sec]')
ylabel('u1')
title('Input u1')
subplot(3,1,3)
stairs(tVec,C*xVec,'LineWidth',2)
hold all;
grid
xlabel('time [sec]')
ylabel('y')
title('Output y')


figure, subplot(3,1,1)
stairs(tVec,[1 0 0 0 0 0]*xVec,'LineWidth',2)
hold all;
grid
xlabel('time [sec]')
ylabel('x')
title('State x')
subplot(3,1,2)
stairs(tVec,[0 1 0 0 0 0]*xVec,'LineWidth',2)
hold all;
grid
xlabel('time [sec]')
ylabel('theta')
title('State theta')
subplot(3,1,3)
stairs(tVec,[0 0 1 0 0 0]*xVec,'LineWidth',2)
hold all;
grid
xlabel('time [sec]')
ylabel('alpha')
title('State alpha')
set(findall(gcf,'Type','line'),'LineWidth',2)
set(findall(gcf,'-property','FontSize'),'FontSize',14);
% legend('$u_{max} = 1.5$','$u_{max} = 2.5$','$u_{max} = 4$')
