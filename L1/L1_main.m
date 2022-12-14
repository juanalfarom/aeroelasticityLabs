%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOMEWORK 1 AEROELASTICITY
% STUDY OF DIVERGENCE
% ISAAC ROBLEDO MARTIN
% FERNANDO RUIZ CERRAJERO
% JUAN ALFARO MORENO
% UC3M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
%% INITIAL DATA

data.L = 10;                      % Structure length [m]
data.E = 72000e06;                % Young modulus [Pa]
data.Cstr = 2.5;                  % Structure chord [m]
data.G = 27100e06;                % Shear modulus [Pa]
data.rho = 2700;                  % Structure density [kg/m^3]
data.nu = 0.1;                    % nu
data.v = 0.33;                    % Poisson`s ratio
data.tskin = @(t) t;              % Upper and lower skin thickness [m]
data.tspar = @(t) 3*t;            % Front and rear spar thickness [m]

%data.tspar = @(t) t;            % Front and rear spar thickness [m]

data.nodes = 8;                   % Number of nodes the beam will be divided in []
data.ntries = 2;                  % Number of tries for the thickness []

%% TASK 01

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumptions: 
%   -Wing is a continuous beam clamped at the wing root
%   -The only contribution to stiffness is from the wing box
%   -The mass is only the structural mass 
% Task:
%   -Estimate t to match results for nu = 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vector with possible thickness [m]s
task1.possible_t = linspace(0.0001,0.003,data.ntries);

task1.firstWB = 3.47;
task1.secondWB = 21.27;
task1.thirdWB = 58.13;
task1.chordwise = 21.75;
task1.torsionWB = 22.57;

task1.firstWB_error = 1e9;
task1.secondWB_error = 1e9;
task1.thirdWB_error = 1e9;
task1.chordwise_error = 1e9;
task1.torsionWB_error = 1e9;
task1.error = 1e9;

% Cross sectional area of the beam as function of thickness [m^2]
task1.A = @(t) (data.nu*data.Cstr^2)-((data.nu*data.Cstr-2*data.tskin(t))*(data.Cstr-2*data.tspar(t)));
% Beam total mass as function of thickness [kg]
task1.Mass = @(t) data.rho*data.L*task1.A(t);
% Second moments of area of the beam as function of thickness [mm^4]
task1.Ixx = @(t) (((1/12)*data.Cstr*(data.nu*data.Cstr)^3)-((1/12)*(data.Cstr-2*data.tspar(t))*(data.nu*data.Cstr-2*data.tskin(t))^3));
task1.Izz = @(t) (((1/12)*data.nu*data.Cstr*(data.Cstr)^3)-((1/12)*(data.nu*data.Cstr-2*data.tskin(t))*(data.Cstr-2*data.tspar(t))^3));
task1.Iyy = @(t) ((1/12)*data.nu*data.Cstr*data.Cstr*(data.Cstr^2 + (data.nu*data.Cstr)^2)) - ((1/12)*(data.nu*data.Cstr-2*data.tskin(t))*(data.Cstr-2*data.tspar(t))*((data.Cstr-2*data.tspar(t))^2 + (data.nu*data.Cstr-2*data.tskin(t))^2));
task1.J   = @(t) 2*3*t^2*(data.Cstr - 3*t)^2*(data.Cstr*data.nu - t)^2/(data.Cstr*3*t + data.nu*data.Cstr*t - 9*t^2 - t^2);
% task1.solver = @(t) (((1/12)*data.nu*data.Cstr*(data.Cstr)^3)-((1/12)*(data.nu*data.Cstr-2*data.tskin(t))*(data.Cstr-2*data.tspar(t))^3))/((data.nu*data.Cstr^2)-((data.nu*data.Cstr-2*data.tskin(t))*(data.Cstr-2*data.tspar(t))));
% eq = @(t)(((task1.firstWB*2*pi*(data.L^2)/((1.875)^2))^2)*data.rho)/data.E - task1.solver(t);

task1.mode1_fun = @(t) ((1.875)^2/(2*pi*data.L^2))*sqrt((data.E/(data.rho*task1.A(t)))*task1.Ixx(t)) - task1.firstWB;
task1.mode2_fun = @(t) ((4.694)^2/(2*pi*data.L^2))*sqrt((data.E/(data.rho*task1.A(t)))*task1.Ixx(t)) - task1.secondWB;
task1.mode3_fun = @(t) ((7.855)^2/(2*pi*data.L^2))*sqrt((data.E/(data.rho*task1.A(t)))*task1.Ixx(t)) - task1.thirdWB;
t0_guess = 0.0001;

task1.t_WB1 = fzero(task1.mode1_fun,t0_guess);
task1.t_WB2 = fzero(task1.mode2_fun,t0_guess);
task1.t_WB3 = fzero(task1.mode3_fun,t0_guess);

% for i = 1:length(possible_t)
% sol = eq(possible_t(i));
% a = abs(sol)
% if a<0.0001
%     break
% end
% end

%If col> dividing_col --> torsion mode
% task1.dividing_col = (2*data.nodes - 2)/2 + 2;
% 
% 
% for i=1:length(task1.possible_t)
%     %Initializing the matrices
%     task1.K = zeros(data.nodes*2,data.nodes*2);
%     task1.masses = zeros(1,data.nodes);
%     task1.M = zeros(data.nodes*2,data.nodes*2);
%     task1.M_minus = zeros(data.nodes*2,data.nodes*2);
% 
%     current_index = 1;
% 
%     for j=1:(data.nodes-1)
%         if j==1||j==(data.nodes-1)
%             l = data.L/(data.nodes-1)/2;
%         else
%             l = data.L/(data.nodes-1);
%         end
% 
%         task1.K(current_index:current_index+3,current_index:current_index+3) = ...
%             task1.K(current_index:current_index+3,current_index:current_index+3) + ...
%                     Stiffness_matrix_beam(data,task1.A(task1.possible_t(i)), ...
%                                           task1.Ixx(task1.possible_t(i)), ...
%                                           task1.Izz(task1.possible_t(i)), ...
%                                           task1.J(task1.possible_t(i)), ...
%                                           data.L/(data.nodes-1));
%         task1.thickness = task1.possible_t(i);
%         task1.area = task1.A(task1.possible_t(i));
%         task1.mass = task1.Mass(task1.possible_t(i));
% 
%         if j==1
%             task1.masses(j) = task1.mass/(data.nodes-1)/2;
%             task1.masses(end) = task1.mass/(data.nodes-1)/2;
%         else
%             task1.masses(j) = task1.mass/(data.nodes-1);
%         end
%         current_index = current_index + 2;
%     end
% 
%     current_index = 1;
% 
%     for j=1:data.nodes
%         task1.M(current_index:current_index+1,current_index:current_index+1) = ...
%                     Mass_matrix_beam(task1.masses(j), ...
%                                      task1.Ixx(task1.possible_t(i)), ...
%                                      task1.Iyy(task1.possible_t(i)), ...
%                                      task1.Izz(task1.possible_t(i)));
%         current_index = current_index + 2;
%     end
% 
%     
%     
%     for j=1:(2*data.nodes)
%         task1.M_minus(j,j) = task1.M(j,j)^(-1/2);
%     end
% 
%     %Clamping
%     task1.M(1,1) = 0;
%     task1.M(2,2) = 0;
%     task1.M_minus(1,1) = 0;
%     task1.M_minus(2,2) = 0;
% 
%     task1.K_changed = task1.M_minus*task1.K*task1.M_minus;
% 
%     [task1.eigenvectors,task1.eigenvalues] = eig(task1.K_changed);
% 
%     task1.freqs = diag(sqrt(task1.eigenvalues))/2/pi;
%     task1.modes = task1.M_minus*(task1.eigenvectors);
% 
%     for k=1:(2*data.nodes)
%         for j=1:data.nodes
%             task1.bending(j,k) = task1.modes((2*j-1),k);
%         end
%         for j=(data.nodes+1):(2*data.nodes)
%             task1.torsion(j-data.nodes,k) = task1.modes(2*(j-data.nodes),k);
%         end
%     end
% 
%     for j=3:task1.dividing_col
%         task1.bending_modes.freqs(j-2) = task1.freqs(j);
%         task1.bending_modes.bending(:,j-2) = task1.bending(:,j);
%         task1.bending_modes.torsion(:,j-2) = task1.torsion(:,j);
%     end
%     
%     for j=(task1.dividing_col+1):(2*data.nodes)
%         task1.torsion_modes.freqs(j-task1.dividing_col) = task1.freqs(j);
%         task1.torsion_modes.bending(:,j-task1.dividing_col) = task1.bending(:,j);
%         task1.torsion_modes.torsion(:,j-task1.dividing_col) = task1.torsion(:,j);
%     end
%     
%     task1.firstWB_error = sqrt(abs(task1.firstWB^2 - task1.bending_modes.freqs(1)^2));
%     task1.secondWB_error = sqrt(abs(task1.secondWB^2 - task1.bending_modes.freqs(2)^2));
%     task1.thirdWB_error = sqrt(abs(task1.thirdWB^2 - task1.bending_modes.freqs(3)^2));
%     task1.chordwise_error = sqrt(abs(task1.chordwise^2 - task1.bending_modes.freqs(1)^2));
%     task1.torsionWB_error = sqrt(abs(task1.torsionWB^2 - task1.torsion_modes.freqs(1)^2));
% 
%     error = task1.firstWB_error + task1.secondWB_error + task1.thirdWB_error + task1.torsionWB_error;
%     
%     if error < task1.error
%         task1.error = error;
%         fprintf('The simulation is closer to converging with t = %f and a total error of %f \n',task1.thickness,task1.error)
%         solution = task1;
%     end
% 
% end
% 
% task1 = solution;
% 
% clear('solution','error','i','j','k')
% 
% figure(1)
% sgtitle(['Bending modes for t = ', num2str(task1.thickness)])
% for i=1:length(task1.bending_modes.freqs)/2+1
%     subplot(data.nodes/2,2,i)
%     plot(task1.bending_modes.bending(:,i))
%     title(['Bending mode with freq = ' num2str(task1.bending_modes.freqs(i))])
% 
%     if i<length(task1.bending_modes.freqs)/2 && data.nodes>10
%         subplot(data.nodes/2,2,i+10)
%         plot(task1.bending_modes.bending(:,i+10))
%         title(['Bending mode with freq = ' num2str(task1.bending_modes.freqs(i+10))])
%     end
% 
% end
% 
% hold off

% figure(2)
% sgtitle(['Torsion modes for t = ', num2str(task1.thickness)])
% for i=1:length(task1.bending_modes.freqs)/2+1
%     subplot(data.nodes/2,2,i)
%     plot(task1.torsion_modes.torsion(:,i))
%     title(['Torsion mode with freq = ' num2str(task1.torsion_modes.freqs(i))])
% 
%     if i<length(task1.bending_modes.freqs)/2 && data.nodes>10
%         subplot(data.nodes/2,2,i+10)
%         plot(task1.torsion_modes.torsion(:,i+10))
%         title(['Torsion mode with freq = ' num2str(task1.torsion_modes.freqs(i+10))])
%     end
% 
% end
% 
% hold off

%% TASK 02

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumptions: 
%   -Wing is a continuous beam clamped at the wing root
%   -The only contribution to stiffness is from the wing box
%   -The mass is only the structural mass 
% Task:
%   -Estimate nu to match results for table 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task1.thickness = 0.0008;

task2.firstWB = 6.74;
task2.secondWB = 40.71;
task2.thirdWB = 108.80;
task2.chordwise = 22.95;
task2.torsion = 38.05;

% Cross sectional area of the beam as function of nu
task2.A = @(nu) (nu*data.Cstr^2)-((nu*data.Cstr-2*data.tskin(task1.thickness))*(data.Cstr-2*data.tspar(task1.thickness)));
% Second moments of area of the beam as function of nu
task2.Ixx = @(nu) (((1/12)*nu*data.Cstr^4)-((1/12)*(data.Cstr-2*data.tspar(task1.thickness))*(nu*data.Cstr-2*data.tskin(task1.thickness))^3));
task2.Izz = @(nu) (((1/12)*nu*data.Cstr^4)-((1/12)*((data.Cstr-2*data.tspar(task1.thickness))^3)*(nu*data.Cstr-2*data.tskin(task1.thickness))));

task2.mode1_fun = @(nu) ((1.875)^2/(2*pi*data.L^2))*sqrt((data.E/(data.rho*task2.A(nu)))*task2.Ixx(nu)) - task2.firstWB;
task2.mode2_fun = @(nu) ((4.694)^2/(2*pi*data.L^2))*sqrt((data.E/(data.rho*task2.A(nu)))*task2.Ixx(nu)) - task2.secondWB;
task2.mode3_fun = @(nu) ((7.855)^2/(2*pi*data.L^2))*sqrt((data.E/(data.rho*task2.A(nu)))*task2.Ixx(nu)) - task2.thirdWB;
%nu0_guess = 0.0001;

nu_vals = linspace(0.1,0.2,5000);
abs_error = 1e36;


% task2.nu_WB1 = fzero(task2.mode1_fun,nu0_guess);
% task2.nu_WB2 = fzero(task2.mode2_fun,nu0_guess);
% task2.nu_WB3 = fzero(task2.mode3_fun,nu0_guess);

for i=1:length(nu_vals)
    error = task2.mode1_fun(nu_vals(i)) + task2.mode2_fun(nu_vals(i)) + task2.mode3_fun(nu_vals(i));
    if error < abs_error
        abs_error = error;
        task2.nu_WB1 = task2.mode1_fun(nu_vals(i));
        task2.nu_WB2 = task2.mode2_fun(nu_vals(i));
        task2.nu_WB3 = task2.mode3_fun(nu_vals(i));
        task2.nu = nu_vals(i);
    end
end
% task2.nu_WB1 = fzero(task2.mode1_fun,nu_guess);
% task2.nu_WB2 = fzero(task2.mode2_fun,nu_guess);
% task2.nu_WB3 = fzero(task2.mode3_fun,nu_guess);

% Call function to solve task 02
%[task2.nu,task2.area] = task2Fcn(data,task2);
% if (task2.nu == 0) || (task2.area <=0)
%     fprintf('No convergence in task 02 \n')
% else
%     fprintf('Convergence in task 02, nu = %f m and cross sectional area = %f m \n',task2.nu,task2.area)
% end

%% TASK 03

task3.t = task1.t_WB3;
task3.J = task1.J(task1.t_WB3);
task3.I = task1.Ixx(task1.t_WB3);
kb = 64*data.E*task3.I/(9*data.L^3);
kt = 4*task3.J*data.G/(3*data.L);


%% TASK 04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumptions: 
%   -Profile chord = 6.25m
%   -FS 15% and RS at 55%. Elastic axis at 35%. Aero force at 25%
%   -Cla = 2pi, e = 0.625m, S = 62.5m^2
% Task:
%   -Obtain incomp. divergence dynamic presure and speed at SL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task4.e = 0.625; 
task4.S = 62.5;
task4.c = 6.25; 
task4.Cl_alpha = 2*pi;
task4.Kt = kt;
task4.EA = 0.35*task4.c;
task4.RS = 0.55*task4.c;
task4.FS = 0.15*task4.c;
task4.AF = 0.25*task4.c;

% ISA data
data.T0 = 288.15; 
data.alpha = -0.0065;
data.rho0 = 1.225;
data.p0 = 101352.9;
data.g = 9.80665;
data.R = 287.1;
data.gamma = 1.4;
data.a = @(h) sqrt(data.gamma*data.R*(data.T0+data.alpha*h));
data.rho = @(h) data.rho0*(1+(data.alpha*h)/data.T0)^(-1-(data.g/(data.R*data.alpha)));
data.T = @(h) data.T0 + data.alpha*h;
data.p = @(h) data.p0*((data.T(h)/data.T0)^(-data.g/(data.alpha*data.R)));

task4.rhoSL = data.rho(0);
task4.aSL = data.a(0);
task4.qdivi = task4.Kt/(task4.S*task4.e*task4.Cl_alpha); 
task4.vSLi = sqrt(task4.qdivi/(0.5*task4.rhoSL)); 
task4.MSL = task4.vSLi/task4.aSL;   

task4.KTAS = task4.vSLi*1.94384; 

fprintf('The incompressible divergence dynamic pressure at sea level qdiv = %f \n', task4.qdivi)
fprintf('The incompressible divergence speed at sea level Vdiv = %f knots \n', task4.KTAS)


%% TASK 05

task5.hRange = linspace(0, 45000, 45000)*0.3048;

for i = 1:length(task5.hRange)
    task5.rho(i) = data.rho(task5.hRange(i));
    task5.a(i) = data.a(task5.hRange(i));
    task5.A = (0.25*(task5.rho(i)^2)*(task5.a(i)^4))/(task4.qdivi^2);
    task5.M(i) = sqrt((-1+sqrt(1+4*task5.A))/(2*task5.A));
    task5.Vkt(i) = task5.M(i)*data.a(task5.hRange(i))*1.94384;
end

task5.KTAS = task5.M(1)*data.a(0)*1.94384;
task4.qdivc = 0.5*data.rho(0)*(task5.M(1)*data.a(0))^2;

fprintf('The compressible divergence dynamic pressure at sea level qdiv = %f \n', task4.qdivc)
fprintf('The compressible divergence speed at sea level Vdiv = %f knots \n', task5.KTAS)

task5.hRange = task5.hRange/0.3048;

figure(1)
plot(task5.M, task5.hRange, 'r--','LineWidth',1.5)
yline(task5.hRange(find(abs(task5.hRange-10000) < 0.3)), 'b','LineWidth',1);
yline(task5.hRange(find(abs(task5.hRange-20000) < 0.5)), 'b','LineWidth',1);
yline(task5.hRange(find(abs(task5.hRange-30000) < 0.5)), 'b','LineWidth',1);
grid on
title('Divergence envelope in M-H plot')
xlabel('Mach number [-]')
ylabel('Height [ft]')
legend('Divergence envelope', 'H1 = 10000ft', 'H2 = 20000ft', 'H3 = 30000ft')

figure(2)
plot(task5.Vkt, task5.hRange, 'r--','LineWidth',1.5)
yline(task5.hRange(find(abs(task5.hRange-10000) < 0.3)), 'b','LineWidth',1);
yline(task5.hRange(find(abs(task5.hRange-20000) < 0.5)), 'b','LineWidth',1);
yline(task5.hRange(find(abs(task5.hRange-30000) < 0.5)), 'b','LineWidth',1);
grid on
title('Divergence envelope in V-H plot')
xlabel('True air speed [kt]')
ylabel('Height [ft]')
legend('Divergence envelope', 'H1 = 10000ft', 'H2 = 20000ft', 'H3 = 30000ft')

%% FUNCTIONS


function [t,area] = task1Fcn(data,task1)
    for i = 1:length(task1.possible_t)
         area = task1.A(task1.possible_t(i));
         C1 = (data.rho*area*data.L*task1.firstWB^2)/(data.E*task1.Ixx(task1.possible_t(i)));
         C2 = (data.rho*area*data.L*task1.secondWB^2)/(data.E*task1.Ixx(task1.possible_t(i)));
         C3 = (data.rho*area*data.L*task1.secondWB^2)/(data.E*task1.Ixx(task1.possible_t(i)));
         
         mode1 = cos(C1)*cosh(C1);
         mode2 = cos(C2)*cosh(C2);
         mode3 = cos(C3)*cosh(C3);
         
         if (mode1==1)&&(mode2==1)&&(mode3==1)
             t = task1.possible_t(i);
             return
         end
    end
    t = 0;
end  % Task 01 solver function

function [nu,area] = task2Fcn(data,task2)
    for i = 1:length(task2.possible_nu)
         area = task2.A(task2.possible_nu(i));
         C1 = (data.rho*area*data.L*task2.firstWB^2)/(data.E*task2.Izz(task2.possible_nu(i)));
         C2 = (data.rho*area*data.L*task2.secondWB^2)/(data.E*task2.Izz(task2.possible_nu(i)));
         C3 = (data.rho*area*data.L*task2.secondWB^2)/(data.E*task2.Izz(task2.possible_nu(i)));
         
         mode1 = cos(C1)*cosh(C1);
         mode2 = cos(C2)*cosh(C2);
         mode3 = cos(C3)*cosh(C3);
         
         if (mode1==1)&&(mode2==1)&&(mode3==1)
             nu = task2.possible_nu(i);
             return
         end
    end
    nu = 0;
end  % Task 02 solver function

