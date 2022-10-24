clear; clc;

%% Task 01-02
t_lin = linspace(0,100,100);
t = linspace(0,100,100);

wagner = @(t) (t+2)./(t+4);
lift_nd = @(t) (t+1)./(t+2);

figure(1)
plot(t_lin,wagner(t_lin),'r',t_lin,lift_nd(t_lin),'b','LineWidth',1.5);
grid on
title('Wagner function \Phi(\tau) vs analytical result L(\tau)')
legend('Wagner function \Phi(\tau)','Analytical result L(\tau)')
xlabel('Semi-chords travelled \tau [-]')
ylabel('Normalised lift [-]')

%% Task 03
factor = 0.5;
figure(2)
plot(t_lin,wagner(t_lin),'r',t_lin,lift_nd(factor*t_lin),'b','LineWidth',1.5);
grid on
title('Wagner function \Phi(\tau) vs analytical result L(\tau)')
legend('Wagner function \Phi(\tau)','Analytical result L(\tau)')
xlabel('Semi-chords travelled \tau [-]')
ylabel('Normalised lift [-]')

%% Task 04
