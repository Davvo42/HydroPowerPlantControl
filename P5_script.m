clc
clear all
close all

T_i=20;
ome_c=0.62394;
T_A=12;
omen=2*50*pi;
Pen=340e6;
BBB=-(T_i/(ome_c*T_A)/(-1*ome_c+T_i))*omen/Pen*-0.1*Pen

S_s = tf([T_A 1 0], [T_i/ome_c T_i+1/ome_c 1])*T_i/ome_c/T_A
dome_dPel = - S_s*tf([omen/Pen], [T_A 1])*(-36e6)

% figure
% bode(S_s)

% figure(1)
% step(dome_dPel)
% hold

% tf([0.1451],[20 1])+tf([-0.1451],[1.6 1]) %%check heavyside

HS_slow = tf([-BBB],[20 1])
HS_fast = tf([BBB],[1.6 1])

figure(2)
step(HS_slow)
hold
step(HS_fast)

figure(1)
step(HS_slow+HS_fast)


% bode(HS_slow)
% hold
% bode(HS_fast)

%% Point 7

Ti = 20;
ome0 = 2*pi*50;
Ta = 0.2;
TA = 12;
T = 0.3507;
kp = 0.006646622;
An = 0.2943;
mu=ome0/(340e6)*-36e6;


Ls = ome0/An*kp * tf([Ti 1], [Ti 0]) * tf([1], [Ta 1]) * tf([-2*T 1], [T 1]) * tf([1], [TA 1])

Ss = 1/(1+Ls)
Ga = -mu*tf([1],[TA 1]) * Ss

figure(2)
margin(Ls)

figure(1)
step(Ga)
hold on;


%figure(1)
%margin(Ls)
%hold

T = 0.3507/5;
TA = 5*12;
kp = 0.006646622;
An = 0.2943/5;
ome0 = 2*pi*50;
%An = 0.29337/5;
mu=ome0/(340e6)*5*-36e6;

Lsnew = ome0/An*kp * tf([Ti 1], [Ti 0]) * tf([1], [Ta 1]) * tf([-2*T 1], [T 1]) * tf([1], [TA 5])

Ssnew = 1/(1+Lsnew)
Ganew = -mu*tf([1],[TA 5]) * Ssnew
step(Ganew)

%margin(Lsnew)
%[Gm, Pm, Wcg, Wcp] = margin(Lsnew)

%% Test