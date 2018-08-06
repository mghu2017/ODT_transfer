% This code is for calculating the atom trajectory in ODT when moving ODT
clear all
close all
clc


d = 300;            %[mm]
f0 = 12;            %[Hz] trapping frequencies
T0 = 1/ f0;
Tperiod = 4*T0;
T1 = Tperiod/2;
T2 = Tperiod - T1;

syms accelx
eqn = 0.5*accelx*T1^2 + accelx*T1*T2 - 0.5*accelx*T2^2 == d;
sol = solve(eqn, accelx);
accel = sol;        %[mm/s^2] acceleration 
decel = accel;      %[mm/s^2] deceleration


xc = @(t) 0.5.*accel.*t.^2.*((t < T1) & (t >= 0)) + (0.5.*accel.*(T1 ...
    ).^2 + accel.*T1.*(t-T1) - 0.5.*decel.*(t-T1).^2).*((t >= T1) & (t <= Tperiod))...
    + 0.*(t < 0) + d.* (t > Tperiod);

t = 0:0.0001:Tperiod;
plot(t, xc(t))

