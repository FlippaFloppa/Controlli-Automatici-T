close all; clear all; clc;

%% Definizioni
syms x1 x2 x3 u N D t opt

b1 = 0.3;
b2 = 0.1;
m=1;
k=1.5;
kG=6.67e-11;
M = 5.98e24;

% Sistema non lineare
f1=x2;
f2=-b1*x2/m + (k-1)*(kG*M/x1^2 -x1*x3^2);
f3=-2*x3*x2/x1 - b2*x3/m + u/(m*x1);
y=x3;

%% Calcolo coppia di equilibrio

x1_e=3e7; % Posizione di equilibrio

F = @(x) [x(2);
         -b1*x(2)/m + (k-1)*(kG*M/x1_e^2 -x1_e*x(3)^2);
         -2*x(3)*x(2)/x1_e - b2*x(3)/m + x(1)/(m*x1_e)];

x=fsolve(F,double([0,0,0]));

x2_e=x(1,2);
x3_e=x(1,3);
u_e=x(1,1);

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

s = tf('s');
[N,D]=ss2tf(A,B,C,D);
G = tf(N,D);

% dimensione grafico
bode_min=1e-10;
bode_max=1e8;

figure(1);
hold on;
grid on;
bode(G,{bode_min,bode_max});
title("Funzione di trasferimento iniziale");
hold off;

%% Definizione parametri regolatore

% ampiezze gradini
WW = 8e-5;
DD = 3e-5;
NN = 2e-4;

% errore a regime
e_star = 0.1;

% attenuazione disturbo sull'uscita
A_d = 45;
omega_d_MAX = 8e-2;
omega_d_min=1e-11;

% attenuazione disturbo di misura
A_n = 85;
omega_n_min = 5e4;
omega_n_MAX = 7.5e7;

% Sovraelongazione massima e tempo d'assestamento al 5%
S_100_spec = 0.01;
T_a5_spec = 0.15;

%% Sintesi Regolatore Statico
mu_s=1e9; % guadagno regolatore statico
Rs=mu_s/(s); % per soddisfare l'errore a regime è stato necessario
             % aggiungere un polo nell'origine
G_e=G*Rs;

%% Grafico le specifiche

figure(2);
hold on;
% Calcolo specifiche S% => Margine di fase
xi = 0.83;  % Coefficiente di smorzamento 0.826085
S_100 = 100*exp(-pi*xi/sqrt(1-xi^2));
Mf_spec = xi*100; % Mf_spec = 83

% Specifiche su d
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -1000; -1000];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 500; 500];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione di taglio)
omega_Ta_min = 1e-11; % lower bound per il plot
omega_Ta_MAX = 460/(Mf_spec*T_a5_spec); % omega_c >= 460/(Mf*T^*) ~ 3.69
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -1000; -1000];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Errore a regime
err_min=1e-11;
err_max=1e-4;
Bnd_err_x=[err_min;err_max;err_max;err_min];
Bnd_err_y=[200;A_d;A_d;A_d];
patch(Bnd_err_x, Bnd_err_y,'m','FaceAlpha',0.2,'EdgeAlpha',0);

% fisica realizzabilità
fis_min=1e6;
fis_max=1e10;
Bnd_fis_x=[fis_min,fis_max,fis_max,fis_max];
Bnd_fis_y=[-A_n,-A_n,-A_n-200,-A_n-200];
patch(Bnd_fis_x, Bnd_fis_y,'y','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}";"e∞";"fis. real.";"G_e(j\omega)"; "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
opt = bodeoptions('cstprefs');
opt.YLim={[-400 400],[-180 0]};
opt.XLim={[bode_min,bode_max]};
bodeplot(G_e,opt);
margin(G,{bode_min,bode_max});
title("Diagramma di Bode con regolatore Statico")
grid on;
zoom on;

% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec- 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_arg = ["G_e(j\omega)";"G(j\omega)"; "M_f"];
legend(Legend_arg);

%% Sintesi Regolatore dinamico

% Rete anticipatrice
Mf_star = Mf_spec; % Mf_star = 83
omega_c_star=1e2;

[mag_omega_c_star, arg_omega_c_star, ~] = bode(G_e, omega_c_star);
mag_omega_c_star_db = 20*log10(mag_omega_c_star);

M_star = 10^(-mag_omega_c_star_db/20);
phi_star = Mf_star - 180 - arg_omega_c_star;

% Formule di inversion
alpha_tau = (cos(phi_star*pi/180) - 1/M_star)/(omega_c_star*sin(phi_star*pi/180))
tau = (M_star - cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180))
Rd_a=(1 + tau*s)/(1 + alpha_tau*s);

check_flag = cos(phi_star*pi/180) - 1/M_star;
if check_flag < 0
    disp('Errore: alpha negativo');
    return;
end

% Regolatore dinamico
Rd = Rd_a;

%% Diagramma di Bode con regolatore Statico e Dinamico

L=Rd*G_e;   % Funzione di anello

figure(3);
hold on;

% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_err_x, Bnd_err_y,'m','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_fis_x, Bnd_fis_y,'y','FaceAlpha',0.2,'EdgeAlpha',0);
Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}";"e∞";"fis. real.";"L(j\omega)"; "G_e(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
bodeplot(L,opt);
margin(G_e,{bode_min,bode_max});
title("Diagramma di Bode con regolatore Statico e Dinamico");
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;
Legend_arg = ["L(j\omega)";"G_e(j\omega)"; "M_f"];
legend(Legend_arg);

%% Definizione funzioni di sensitività
R=Rd*Rs;    % Regolatore

S=1/(1+R*G); %  Funzione di sensitività
F=(R*G)/(1+R*G); %  Funzione di sensitività complementare
Q=R/(1+R*G); %  Funzione di sensitività del controllo

figure(10)
hold on;grid on;
legend(["S(j\omega)";"F(j\omega)"; "Q(j\omega)"]);
margin(S,{bode_min,bode_max});
margin(F,{bode_min,bode_max});
margin(Q,{bode_min,bode_max});
legend(["S(j\omega)";"F(j\omega)"; "Q(j\omega)"]);
title("Funzioni di sensitività relative a L(j\omega)");

%% Risposta al gradino
figure(4);
% disturbo di ingresso
W=8e-5; 
T_simulation = 1;
[y_step,t_step] = step(W*F, 1);
grid on, zoom on, hold on;

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[W*(1+S_100_spec),W*(1+S_100_spec),W+1,W+1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
patch([T_a5_spec,T_simulation,T_simulation,T_a5_spec],[W-W*0.05,W-W*0.05,0,0],'g','FaceAlpha',0.3,'EdgeAlpha',0.5)
ylim([W/2,W+1]);

plot(t_step,y_step,'b');
title("Risposta al gradino");
ylim([0,W*1.2]);

Legend_step = ["Vincolo sovraelongazione";"Tempo di assestamento";"Risposta al gradino"];
legend(Legend_step);

%% Check disturbo in uscita

figure(5);

%intervallo di campionamento
tt = (0:1e-2:1e3);

% calcolo disturbo di uscita
d=0;
for i=1:4
    d=d+sin(0.02*i*tt);
end
d=DD*d;

%plot
y_n = lsim(S,d,tt);
hold on, grid on, zoom on
plot(tt,d,'m')
plot(tt,y_n,'b')
title("Comportamento disturbo di uscita");
grid on
legend('d','y_d')

%% Check disturbo di misura
figure(6);

%intervallo di campionamento
tt = (0:1e-6:1e-2);
% calcolo disturbo di misura
n=0;
for i=1:4
    n=n+sin(5e4*i*tt);
end
n=NN*n;

% plot
y_n = lsim(-F,n,tt);
hold on, grid on, zoom on
plot(tt,n,'m');
plot(tt,y_n,'b');
title("Comportamento disturbo di misura");
grid on
legend('n','y_n');