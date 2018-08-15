% This code is for calculating the atom trajectory in ODT when moving ODT
clear all
clc

d = 315/2.727;            %   [mm] tansfer distance
f0 = 10;            %   [Hz] trapping frequencies
T0 = 1/ f0;
w0 = 2*pi*f0;

%%%%%%--------Numerical way------------
% define analytic atom trajectory in lab frame
syms accAvg(dist, Ttrans)
accAvg(dist, Ttrans) = 4*dist/Ttrans^2;    % [mm/s^2]
% define acceleration of ODT transfer
syms accelFunc(tx, Ttrans)
accelFunc(tx, Ttrans) = (2*accAvg(d, Ttrans)).*(triangularPulse(0,Ttrans/4,Ttrans/2,tx)-triangularPulse(Ttrans/2,3*Ttrans/4,Ttrans,tx));
% define velocity profile of ODT transfer
syms velFunc(tx, Ttrans) tx1
velFunc(tx, Ttrans) = int(accelFunc(tx1, Ttrans), tx1, 0, tx);
% define ODT trajectory
syms disFunc(tx, Ttrans)
disFunc(tx, Ttrans) = int(velFunc(tx, Ttrans), tx, 0, tx);
% define atom trajectory in lab frame
syms x(tx, Ttrans) t
x(tx, Ttrans) = disFunc(tx, Ttrans) + (1/2/pi/f0).*int(sin(2.*pi.*f0.*(t-tx)).*accelFunc(t, Ttrans), t,...
    0, tx);
% define atom slosh amplitude
syms A(Ttrans, w) tx
A(Ttrans, w) = abs(fourier(velFunc(tx, Ttrans), tx, w));


close all
t2 = 8*T0;
t1 = (-0.1:0.01:(t2/T0)+0.1).*T0;
h1 = figure();
Nline = 3;
Ncol = 3;
set(h1, 'Position', [0 100 800 800]);         %[left bottom width height]
subplot(Nline,Ncol,1);
plot(t1./T0,accelFunc(t1, t2));
xlabel('t (T0)');
ylabel('Acceleration (mm/s^2)');
accAvg0 = double(accAvg(d, t2));                    % [mm/s^2]
title(['Avg accel = ', num2str(accAvg0), ' mm/s^2 (Trans = ', num2str(t2/T0), 'T0)']);
subplot(Nline,Ncol,2);
vel = double(velFunc(t1, t2));
plot(t1./T0, vel);
xlabel('t (T0)');
ylabel('Velocity (mm/s)');
title(['peak velocity = ', num2str(max(vel)), ' mm/s']);
subplot(Nline,Ncol,3);
plot(t1./T0, disFunc(t1, t2));
xlabel('t (T0)');
ylabel('x_c(t) (mm)');
title(['trap period T_0 = ', num2str(T0*1000), ' ms']);
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
title(['Atom slosh in lab frame (Ttrans = ', num2str(t2/T0),'T0)']);
subplot(Nline, Ncol, 8);
plot(t1./T0, x(t1, t2)-disFunc(t1,t2));
hold on
plot(t3./T0, x(t3, t2)-disFunc(t3,t2), 'r', 'LineWidth', 2);
grid on
grid minor
xlabel('t (T0)');
ylabel('x(t)-x_c(t) (mm)');
title(['Atom slosh in trap frame (Ttrans = ', num2str(t2/T0),'T0)']);
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
title(['Atom slosh in trap frame (Ttrans = ', num2str(t2/T0),'T0)']);


subplot(Nline,Ncol,4);
Ttransx = 5*T0:0.02:(15*T0);
plot(Ttransx./T0, A(Ttransx, w0));
xlabel('Ttrans (T0)');
ylabel('Slosh amplitude (mm)');

subplot(Nline,Ncol,5);
t2 = 3*T0;
t1 = t2:0.002:(t2+3*T0);
plot(t1./T0, x(t1, t2));
xlabel('t (T0)');
ylabel('Slosh (mm)');
title(['Slosh at Ttrans = ', num2str(t2/T0), 'T0']);

subplot(Nline,Ncol,6);
t2 = 4*T0;
t1 = t2:0.002:(t2+3*T0);
plot(t1./T0, x(t1, t2));
xlabel('t (T0)');
ylabel('Slosh (mm)');
title(['Slosh at Ttrans =', num2str(t2/T0), 'T0']);

