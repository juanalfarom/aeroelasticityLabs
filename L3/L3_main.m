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

test_velocities = linspace(215,275,150);

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

flutter.glo_linear = 100;
flutter.idxlo_linear = 0;
for i=1:length(test_velocities)
    if (abs(mode2.glo_linear(i))<flutter.glo_linear)&&(mode2.glo_linear(i)>0)
        flutter.glo_linear = mode2.g_linear(i);
        flutter.idxlo_linear = i;
    end
end

flutter.Vlo_linear = test_velocities(flutter.idxlo_linear);

flutter.gup_linear = 100;
flutter.idxup_linear = 0;
for i=1:length(test_velocities)
    if (abs(mode2.gup_linear(i))<flutter.gup_linear)&&(mode2.gup_linear(i)>0)
        flutter.gup_linear = mode2.g_linear(i);
        flutter.idxup_linear = i;
    end
end

flutter.Vup_linear = test_velocities(flutter.idxup_linear);

%% Section 4

mode1.w = mode1.f*2*pi;
mode2.w = mode2.f*2*pi;

mode1.wlo = mode1.flo*2*pi;
mode2.wlo = mode2.flo*2*pi;

mode1.wup = mode1.fup*2*pi;
mode2.wup = mode2.fup*2*pi;

mode1.beta = mode1.w.*mode1.g;
mode2.beta = mode2.w.*mode2.g;

mode1.betalo = mode1.wlo.*mode1.glo;
mode2.betalo = mode2.wlo.*mode2.glo;

mode1.betaup = mode1.wup.*mode1.gup;
mode2.betaup = mode2.wup.*mode2.gup;


flutter.F = ((mode2.w.^2 - mode1.w.^2)/2 + (mode2.beta.^2 - mode1.beta.^2)/2).^2 + ...
            4.*mode1.beta.*mode2.beta.*((mode2.w.^2 + mode1.w.^2)/2 + 2*(mode2.beta + mode1.beta).^2/4) - ...
            (((mode2.beta - mode1.beta)./(mode2.beta + mode1.beta)).*(mode2.w.^2 - mode1.w.^2)/2 + 2*(mode2.beta + mode1.beta).^2/4).^2;

flutter.Flo = ((mode2.wlo.^2 - mode1.wlo.^2)/2 + (mode2.betalo.^2 - mode1.betalo.^2)/2).^2 + ...
            4.*mode1.betalo.*mode2.betalo.*((mode2.wlo.^2 + mode1.wlo.^2)/2 + 2*(mode2.betalo + mode1.betalo).^2/4) - ...
            (((mode2.betalo - mode1.betalo)./(mode2.betalo + mode1.betalo)).*(mode2.wlo.^2 - mode1.wlo.^2)/2 + 2*(mode2.betalo + mode1.betalo).^2/4).^2;

flutter.Fup = ((mode2.wup.^2 - mode1.wup.^2)/2 + (mode2.betaup.^2 - mode1.betaup.^2)/2).^2 + ...
            4.*mode1.betaup.*mode2.betaup.*((mode2.wup.^2 + mode1.wup.^2)/2 + 2*(mode2.betaup + mode1.betaup).^2/4) - ...
            (((mode2.betaup - mode1.betaup)./(mode2.betaup + mode1.betaup)).*(mode2.wup.^2 - mode1.wup.^2)/2 + 2*(mode2.betaup + mode1.betaup).^2/4).^2;

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

%

flutter.fitlo = polyfit(q,flutter.Flo,2);
flutter.Flo_approx = polyval(flutter.fitlo,test_q);

flutter.Flo_flutter = 100;
flutter.Flo_idx = 1;

for i=1:length(flutter.Flo_approx)
    if (abs(flutter.Flo_approx(i))<flutter.Flo_flutter)&&(flutter.Flo_approx(i)>0)
        flutter.Flo_flutter = flutter.Flo_approx(i);
        flutter.Flo_idx = i;
    end
end

flutter.Vlo_flutter = sqrt(test_q(flutter.Flo_idx)/(0.5*1.225));
flutter.qlo_flutter = test_q(flutter.Flo_idx);

%

flutter.fitup = polyfit(q,flutter.Fup,2);
flutter.Fup_approx = polyval(flutter.fitup,test_q);

flutter.Fup_flutter = 100;
flutter.Fup_idx = 1;

for i=1:length(flutter.Fup_approx)
    if (abs(flutter.Fup_approx(i))<flutter.Fup_flutter)&&(flutter.Fup_approx(i)>0)
        flutter.Fup_flutter = flutter.Fup_approx(i);
        flutter.Fup_idx = i;
    end
end

flutter.Vup_flutter = sqrt(test_q(flutter.Fup_idx)/(0.5*1.225));
flutter.qup_flutter = test_q(flutter.Fup_idx);

%

flutter.Ferror = mean((abs(flutter.F-flutter.Flo) + abs(flutter.F - flutter.Fup))/2./flutter.F);
flutter.F_lo_error = (flutter.F_approx(1)*(1-flutter.Ferror) - flutter.F_approx(1)) + flutter.F_approx;
flutter.F_lo_flutter = 100;
for i=1:length(flutter.F_lo_error)
    if (abs(flutter.F_lo_error(i))<flutter.F_lo_flutter)&&(flutter.F_lo_error(i)>0)
        flutter.F_lo_flutter = flutter.F_lo_error(i);
        flutter.F_lo_idx = i;
    end
end
flutter.V_low_flutter = sqrt(test_q(flutter.F_lo_idx)/(0.5*1.225));
flutter.q_low_flutter = test_q(flutter.F_lo_idx);

flutter.F_up_error = (flutter.F_approx(1)*(1+flutter.Ferror) - flutter.F_approx(1)) + flutter.F_approx;

flutter.F_up_flutter = 100;
for i=1:length(flutter.F_up_error)
    if (abs(flutter.F_up_error(i))<flutter.F_up_flutter)&&(flutter.F_up_error(i)>0)
        flutter.F_up_flutter = flutter.F_up_error(i);
        flutter.F_up_idx = i;
    end
end
flutter.V_up_flutter = sqrt(test_q(flutter.F_up_idx)/(0.5*1.225));
flutter.q_up_flutter = test_q(flutter.F_up_idx);

%% Figures
set(groot, 'defaultLegendFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultAxesLineWidth', 1);
set(groot, 'defaultAxesXMinorTick', 'on');
set(groot, 'defaultAxesYMinorTick', 'on');
set(groot, 'defaultLegendBox', 'off');
set(groot, 'defaultLegendLocation', 'best');
set(groot, 'defaultLineLineWidth', 1);
set(groot, 'defaultLineMarkerSize', 5);
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

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
    subplot(1,2,1)
    plot(velocities,flutter.F,'b',test_velocities,flutter.F_approx,'b--')
    hold on
    plot(velocities,flutter.Flo,'r',test_velocities,flutter.Flo_approx,'r--')
    hold on
    plot(velocities,flutter.Fup,'r',test_velocities,flutter.Fup_approx,'r--')
    hold on
    yline(0,'k--')
    xline(flutter.V_flutter,'k--',{'Flutter KCAS'});
    ylabel('F [-]','fontsize',14,'interpreter','latex')
    xlabel('Airspeed [Kts]','fontsize',14,'interpreter','latex')
    legend('Data','Approximation','Location','southwest','fontsize',14,'interpreter','latex')
    title('Zimmerman & Weissenburger approximation','fontsize',14,'interpreter','latex')

    subplot(1,2,2)
    plot(velocities,flutter.F,'b',test_velocities,flutter.F_approx,'b--')
    hold on
    plot(velocities,flutter.Flo,'r',test_velocities,flutter.F_lo_error,'r--')
    hold on
    plot(velocities,flutter.Fup,'r',test_velocities,flutter.F_up_error,'r--')
    hold on
    yline(0,'k--')
    xline(flutter.V_flutter,'k--',{'Flutter KCAS'});
    ylabel('F [-]','fontsize',14,'interpreter','latex')
    xlabel('Airspeed [m/s]','fontsize',14,'interpreter','latex')
    legend('Data','Approximation','Location','southwest','fontsize',14,'interpreter','latex')
    title('Zimmerman & Weissenburger approximation','fontsize',14,'interpreter','latex')
end