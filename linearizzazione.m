close all; clear all; clc;

%% Definizioni
syms x1 x2 x3 u b1 b2 m k kG M autovettori autovalori T_inv A_hat B_hat C_hat D_hat modello

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

A = jacobian([f1,f2,f3],[x1,x2,x3])
B = jacobian([f1,f2,f3],u)
C = jacobian(y,[x1,x2,x3])
D = jacobian(y,u)

[T_inv,A_hat]=jordan(A);
A_hat=double(subs(A_hat,[x1,x2,x3,u],[x_e,u_e]))
B_hat = double(subs(inv(T_inv)*B,[x1,x2,x3,u],[x_e,u_e]))
C_hat = double(subs(C*T_inv,[x1,x2,x3,u],[x_e,u_e]))
D_hat = double(subs(D,[x1,x2,x3,u],[x_e,u_e]))

%% Trovo sistema lineare

x0=[x1_e;x2_e;x3_e];
modello = ss(A_hat, B_hat, C_hat, D_hat);
x0_J = inv(T_inv)*x0;   % stato iniziale nelle nuove coordinate
x0_J= double(subs(x0_J,[x1,x2,x3,u],[x_e,u_e]))

uu = zeros(length(interv), 1); % input nullo (evoluzione libera)
[YY_J, TT_J, XX_J] = lsim(modello, uu, interv, x0_J);
YY_J = real(YY_J); % convertiamo il tipo di dato in 'reale' (la parte immaginaria è sicuramente nulla)

%% Grafico
figure;
plot(TT_J,YY_J)
hold on; grid on; zoom on; box on;
title('Traiettoria sistema')
xlim([0, 100])
xlabel('tempo [s]')
ylabel('posizione')
