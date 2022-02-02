close all; clear all; clc;

%% Definizioni
syms x1 x2 x3 u N D

b1 = 0.3;
b2 = 0.1;
m=1;
k=1.5;
kG=6.67e-11;
M = 5.98e24;

% intervallo di tempo
interv = 0:0.1:100; % da 0 a 100 secondi con passo 0.1

% sistema
f1=x2;
f2=-b1*x2/m + (k-1)*(kG*M/x1^2 -x1*x3^2);
f3=-2*x3*x2/x1 - b2*x3/m + u/(m*x1);
y=x3;

%% Calcolo coppia di equilibrio

x1_e=3e7 % Posizione di equilibrio

F = @(x) [x(2);
         -b1*x(2)/m + (k-1)*(kG*M/x1_e^2 -x1_e*x(3)^2);
         -2*x(3)*x(2)/x1_e - b2*x(3)/m + x(1)/(m*x1_e)];

x=fsolve(F,double([0,0,0]));

x2_e=x(1,2)
x3_e=x(1,3)
u_e=x(1,1)

%{
% Punti di equilibrio

x1_e=3e7;
x2_e=0;
x3_e=0.000121543;
u_e=364.63;

%}

% x_e,u_e coppia di equilibrio
x_e=[x1_e,x2_e,x3_e];

%% Jacobiana e linearizzazione

A = jacobian([f1,f2,f3],[x1,x2,x3]);
B = jacobian([f1,f2,f3],u);
C = jacobian(y,[x1,x2,x3]);
D = jacobian(y,u);

A = double(subs(A,[x1,x2,x3,u],[x_e,u_e]));
B = double(subs(B,[x1,x2,x3,u],[x_e,u_e]));
C = double(subs(C,[x1,x2,x3,u],[x_e,u_e]));
D = double(subs(D,[x1,x2,x3,u],[x_e,u_e]));

%% Trasformata di Laplace

figure;
s = tf('s');
[N,D]=ss2tf(A,B,C,D);
G = tf(N,D);
bode(G);
zpk(G)
