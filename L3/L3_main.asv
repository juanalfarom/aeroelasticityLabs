clc;clear;close all;



show_figures = 1;

g1 = 0.008;
g2 = 0.036;

velocities = [40,105,180,215];

mode1.f = [5.42,5.4,5.4,5.075];
mode1.scatter_f = [0.01,0.02,0.03,0.04];
mode2.f = [4.6,4.65,4.95,4.65];
mode2.scatter_f = [0.01,0.02,0.03,0.04];

mode1.g = [g1+0.005,g1+0.01,g1+0.015,g1+0.02];
mode1.scatter_g = [0.02,0.04,0.06,0.08];
mode2.g = [g2/4,g2/2,g2,g2/2];
mode2.scatter_g = [0.02,0.04,0.06,0.08];

%% Section 2

mode1.fup = mode1.f + mode1.scatter_f.*mode1.f;
mode1.flo = mode1.f - mode1.scatter_f.*mode1.f;

mode2.fup = mode2.f + mode2.scatter_f.*mode2.f;
mode2.flo = mode2.f - mode2.scatter_f.*mode2.f;

mode1.gup = mode1.g + mode1.scatter_g.*mode1.g;
mode1.glo = mode1.g - mode1.scatter_g.*mode1.g;

mode2.gup = mode2.g + mode2.scatter_g.*mode2.g;
mode2.glo = mode2.g - mode2.scatter_g.*mode2.g;

test_velocities = linspace(220,275,150);

mode1.g_linear = spline(velocities(2:end),mode1.g(2:end),test_velocities);
mode1.gup_linear = spline(velocities(2:end),mode1.gup(2:end),test_velocities);
mode1.glo_linear = spline(velocities(2:end),mode1.glo(2:end),test_velocities);

mode1.f_linear = spline(velocities(2:end),mode1.f(2:end),test_velocities);
mode1.fup_linear = spline(velocities(2:end),mode1.fup(2:end),test_velocities);
mode1.flo_linear = spline(velocities(2:end),mode1.flo(2:end),test_velocities);

mode2.g_linear = spline(velocities(2:end),mode2.g(2:end),test_velocities);
mode2.gup_linear = spline(velocities(2:end),mode2.gup(2:end),test_velocities);
mode2.glo_linear = spline(velocities(2:end),mode2.glo(2:end),test_velocities);

mode2.f_linear = spline(velocities(2:end),mode2.f(2:end),test_velocities);
mode2.fup_linear = spline(velocities(2:end),mode2.fup(2:end),test_velocities);
mode2.flo_linear = spline(velocities(2:end),mode2.flo(2:end),test_velocities);

flutter.g_linear = 100;
flutter.idx_linear = 0;

for i=1:length(test_velocities)
    if (abs(mode2.g_linear(i))<flutter.g_linear)&&(mode2.g_linear(i)>0)
        flutter.g_linear = mode2.g_linear(i);
        flutter.idx_linear = i;
    end
end

flutter.V_linear = test_velocities(flutter.idx_linear);
flutter.glo_linear = mode2.glo_linear(flutter.idx_linear);

%% Section 4

mode1.w = mode1.f*2*pi;
mode2.w = mode2.f*2*pi;

mode1.beta = mode1.w.*mode1.g;
mode2.beta = mode2.w.*mode2.g;


flutter.F = ((mode2.w.^2 - mode1.w.^2)/2 + (mode2.beta.^2 - mode1.beta.^2)/2).^2 + ...
            4.*mode1.beta.*mode2.beta.*((mode2.w.^2 + mode1.w.^2)/2 + 2*(mode2.beta + mode1.beta).^2/4) - ...
            (((mode2.beta - mode1.beta)./(mode2.beta + mode1.beta)).*(mode2.w.^2 - mode1.w.^2)/2 + 2*(mode2.beta + mode1.beta).^2/4).^2;

q = 0.5*1.225.*velocities.^2;
test_q = 0.5*1.225.*test_velocities.^2;

flutter.fit = polyfit(q,flutter.F,2);
flutter.F_approx = polyval(flutter.fit,test_q);

flutter.F_flutter = 100;
flutter.F_idx = 1;

for i=1:length(flutter.F_approx)
    if (abs(flutter.F_approx(i))<flutter.F_flutter)&&(flutter.F_approx(i)>0)
        flutter.F_flutter = flutter.F_approx(i);
        flutter.F_idx = i;
    end
end

flutter.V_flutter = sqrt(test_q(flutter.F_idx)/(0.5*1.225));
flutter.q_flutter = test_q(flutter.F_idx);

%% figures

if show_figures == 1
    figure (1)
    subplot(2,2,1)
    errorbar(velocities,mode1.f,mode1.scatter_f.*mode1.f,'-bo','MarkerFaceColor','b')
    hold on
    errorbar(velocities,mode2.f,mode2.scatter_f.*mode2.f,'-ro','MarkerFaceColor','r')
    hold on
    plot(velocities,mode1.fup,'b--',velocities,mode1.flo,'b--',velocities,mode2.fup,'r--',velocities,mode2.flo,'r--')
    ylim([3,7])
    xlim([0,240])
    grid on
    xlabel('V (KCAS)')
    ylabel('Frequency [Hz]')
    legend('Mode 1','Mode 2','Location','best')
    title('Frequency Data')

    subplot(2,2,2)
    errorbar(velocities,mode1.g,mode1.scatter_g.*mode1.g,'-bo','MarkerFaceColor','b')
    hold on
    errorbar(velocities,mode2.g,mode2.scatter_g.*mode2.g,'-ro','MarkerFaceColor','r')
    hold on
    plot(velocities,mode1.gup,'b--',velocities,mode1.glo,'b--',velocities,mode2.gup,'r--',velocities,mode2.glo,'r--')
    ylim([-0.003,max(mode2.g)+0.01])
    xlim([0,240])
    yline(0,'k--')
    grid on
    xlabel('V (KCAS)')
    ylabel('Damping [-]')
    legend('Mode 1','Mode 2','Location','best')
    title('Damping Data')

    subplot(2,2,3)
    plot([velocities,test_velocities],[mode1.f,mode1.f_linear],'-b')
    hold on
    plot([velocities,test_velocities],[mode2.f,mode2.f_linear],'-r')
    hold on
    plot([velocities,test_velocities],[mode1.fup,mode1.fup_linear],'b--',[velocities,test_velocities],[mode1.flo,mode1.flo_linear],'b--',[velocities,test_velocities],[mode2.fup,mode2.fup_linear],'r--',[velocities,test_velocities],[mode2.flo,mode2.flo_linear],'r--')
    ylim([3,7])
    xlim([0,max(test_velocities)+10])
    grid on
    xlabel('V (KCAS)')
    ylabel('Frequency [Hz]')
    xline(flutter.V_linear,'k--',{'Flutter speed'});
    legend('Mode 1','Mode 2','Location','best')
    title('Direct Frequency Calculations')

    subplot(2,2,4)
    plot([velocities,test_velocities],[mode1.g,mode1.g_linear],'-b')
    hold on
    plot([velocities,test_velocities],[mode2.g,mode2.g_linear],'-r')
    hold on
    plot([velocities,test_velocities],[mode1.gup,mode1.gup_linear],'b--',[velocities,test_velocities],[mode1.glo,mode1.glo_linear],'b--',[velocities,test_velocities],[mode2.gup,mode2.gup_linear],'r--',[velocities,test_velocities],[mode2.glo,mode2.glo_linear],'r--')
    hold on
    plot(flutter.V_linear,flutter.g_linear,'k*')
    xlim([0,max(test_velocities)+10])
    yline(0,'k--')
    xline(flutter.V_linear,'k--',{'Flutter speed'});
    grid on
    xlabel('V (KCAS)')
    ylabel('Damping [-]')
    legend('Mode 1','Mode 2','Location','best')
    title('Direct Damping Calculations')

    figure(2)
    plot(q,flutter.F,'b',test_q,flutter.F_approx,'b--')
    yline(0,'k--')
    xline(flutter.q_flutter,'k--',{'Flutter dynamic pressure'});
    ylabel('F [-]')
    xlabel('Dynamic pressure [kg/ms^2]')
    legend('Data','Approximation','Location','southwest')
    title('Zimmerman & Weissenburger approximation')
end