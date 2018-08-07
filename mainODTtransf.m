% This code is for calculating the atom trajectory in ODT when moving ODT
clear all
close all
clc


d = 300;            %[mm]
f0 = 12;            %[Hz] trapping frequencies
T0 = 1/ f0;
Tperiod = 3.5*T0;
T1 = Tperiod/2;
T2 = Tperiod - T1;

syms accelx
eqn = 0.5*accelx*T1^2 + accelx*T1*T2 - 0.5*accelx*T2^2 == d;
sol = solve(eqn, accelx);
accel = sol;        %[mm/s^2] acceleration 
decel = accel;      %[mm/s^2] deceleration

%%%%%%--------Numerical way------------
syms tx
accelFunc(tx) = piecewise(0<=tx<=T1, accel, T1<tx<=Tperiod, - decel, ...
    tx<0, 0, tx>Tperiod, 0);
syms velFunc(tx)
velFunc(tx) = int(accelFunc(tx), tx, 0, tx);
syms disFunc(tx)
disFunc(tx) = int(velFunc(tx), tx, 0, tx);
syms x(tx) t
x(tx) = disFunc(tx) + (1/2/pi/f0).*int(sin(2.*pi.*f0.*(t-tx)).*accelFunc(t), t,...
    0, tx);

t1 = 0:0.002:Tperiod;
h1 = figure();
set(h1, 'Position', [50 100 500 800]);         %[left bottom width height]
subplot(4,1,1);
plot(t1,accelFunc(t1));
subplot(4,1,2);
plot(t1, velFunc(t1));
subplot(4,1,3);
plot(t1, disFunc(t1));
subplot(4,1,4);
t2 = Tperiod:0.002:3*Tperiod;
plot(t2, x(t2));
