clc;clear;


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

mode1.fup = mode1.f + mode1.scatter_f.*mode1.f;
mode1.flo = mode1.f - mode1.scatter_f.*mode1.f;

mode2.fup = mode2.f + mode2.scatter_f.*mode2.f;
mode2.flo = mode2.f - mode2.scatter_f.*mode2.f;

mode1.gup = mode1.g + mode1.scatter_g.*mode1.g;
mode1.glo = mode1.g - mode1.scatter_g.*mode1.g;

mode2.gup = mode2.g + mode2.scatter_g.*mode2.g;
mode2.glo = mode2.g - mode2.scatter_g.*mode2.g;



if show_figures == 1
    figure (1)
    subplot(1,2,1)
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
    legend('Mode 1','Mode 2')
    title('Frequency Data')

    subplot(1,2,2)
    errorbar(velocities,mode1.g,mode1.scatter_g.*mode1.g,'-bo','MarkerFaceColor','b')
    hold on
    errorbar(velocities,mode2.g,mode2.scatter_g.*mode2.g,'-ro','MarkerFaceColor','r')
    hold on
    plot(velocities,mode1.gup,'b--',velocities,mode1.glo,'b--',velocities,mode2.gup,'r--',velocities,mode2.glo,'r--')
    ylim([-0.003,max(mode2.g)+0.01])
    xlim([0,240])
    grid on
    xlabel('V (KCAS)')
    ylabel('Damping [-]')
    legend('Mode 1','Mode 2')
    title('Damping Data')
end