% This code is for calculating the atom trajectory in ODT when moving ODT
clear all
clc

d = 300;            %   [mm] tansfer distance
f0 = 12;            %   [Hz] trapping frequencies
T0 = 1/ f0;
w0 = 2*pi*f0;

syms accel(Tperiod) decel(Tperiod)
accel(Tperiod) = d/(Tperiod/2)^2;     %[mm/s^2] acceleration 
decel(Tperiod) = d/(Tperiod/2)^2;     %[mm/s^2] deceleration

%%%%%%--------Numerical way------------
syms accelFunc(tx, Tperiod)
accelFunc(tx, Tperiod) = piecewise(0<=tx<=(Tperiod/2), accel(Tperiod), (Tperiod/2)<tx<=Tperiod, - decel(Tperiod), ...
    tx<0, 0, tx>Tperiod, 0);
% syms velFunc(tx, Tperiod)
% velFunc(tx, Tperiod) = int(accelFunc(tx, Tperiod), tx, 0, tx);
syms velFunc(tx, Tperiod)
velFunc(tx, Tperiod) = accel(Tperiod).*(Tperiod/2).*triangularPulse(0,Tperiod/2,Tperiod,tx);
syms disFunc(tx, Tperiod)
disFunc(tx, Tperiod) = int(velFunc(tx, Tperiod), tx, 0, tx);
syms x(tx, Tperiod) t
x(tx, Tperiod) = disFunc(tx, Tperiod) + (1/2/pi/f0).*int(sin(2.*pi.*f0.*(t-tx)).*accelFunc(t, Tperiod), t,...
    0, tx);
syms A(Tperiod, w) tx
A(Tperiod, w) = abs(fourier(velFunc(tx, Tperiod), tx, w));

close all
t2 = 4*T0;
t1 = (-0.1:0.01:(t2/T0)+0.1).*T0;
h1 = figure();
Nline = 3;
Ncol = 3;
set(h1, 'Position', [0 100 800 800]);         %[left bottom width height]
subplot(Nline,Ncol,1);
plot(t1./T0,accelFunc(t1, t2));
xlabel('t (T0)');
ylabel('Acceleration (mm/s^2)');
subplot(Nline,Ncol,2);
plot(t1./T0, velFunc(t1, t2));
xlabel('t (T0)');
ylabel('Velocity (mm/s)');
subplot(Nline,Ncol,3);
plot(t1./T0, disFunc(t1, t2));
xlabel('t (T0)');
ylabel('x_c(t) (mm)');

subplot(Nline, Ncol, 7);
t2 = 3.*T0;
t1 = (-0.1:0.01:(t2/T0)).*T0;
t3 = (t2/T0:0.01:(t2/T0)+2).*T0;
plot(t1./T0, x(t1, t2));
hold on
plot(t3./T0, x(t3, t2), 'r', 'LineWidth', 2);
grid on
grid minor
xlabel('t (T0)');
ylabel('x(t) (mm)');
title(['Atom slosh in lab frame (Tperiod = ', num2str(t2/T0),'T0)']);
subplot(Nline, Ncol, 8);
plot(t1./T0, x(t1, t2)-disFunc(t1,t2));
hold on
plot(t3./T0, x(t3, t2)-disFunc(t3,t2), 'r', 'LineWidth', 2);
grid on
grid minor
xlabel('t (T0)');
ylabel('x(t)-x_c(t) (mm)');
title(['Atom slosh in trap frame (Tperiod = ', num2str(t2/T0),'T0)']);
subplot(Nline, Ncol, 9);
t2 = 4*T0;
t1 = (-0.1:0.01:(t2/T0)).*T0;
t3 = (t2/T0:0.01:(t2/T0)+2).*T0;
plot(t1./T0, x(t1, t2)-disFunc(t1,t2));
hold on
plot(t3./T0, x(t3, t2)-disFunc(t3,t2), 'r', 'LineWidth', 2);
grid on
grid minor
xlabel('t (T0)');
ylabel('x(t)-x_c(t) (mm)');
title(['Atom slosh in trap frame (Tperiod = ', num2str(t2/T0),'T0)']);


subplot(Nline,Ncol,4);
Tperiodx = 0.01*T0:0.002:(5*T0);
plot(Tperiodx./T0, A(Tperiodx, w0));
xlabel('Tperiod (T0)');
ylabel('Slosh amplitude (mm)');

subplot(Nline,Ncol,5);
t2 = 3*T0;
t1 = t2:0.002:(t2+3*T0);
plot(t1./T0, x(t1, t2));
xlabel('t (T0)');
ylabel('Slosh (mm)');
title(['Slosh at Tperiod = ', num2str(t2/T0), 'T0']);

subplot(Nline,Ncol,6);
t2 = 4*T0;
t1 = t2:0.002:(t2+3*T0);
plot(t1./T0, x(t1, t2));
xlabel('t (T0)');
ylabel('Slosh (mm)');
title(['Slosh at Tperiod =', num2str(t2/T0), 'T0']);

