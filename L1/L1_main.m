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
data.nodes = 8;                   % Number of nodes the beam will be divided in []
data.ntries = 1;                  % Number of tries for the thickness []

%% TASK 01

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumptions: 
%   -Wing is a continuous beam clamped at the wing root
%   -The only contribution to stiffness is from the wing box
%   -The mass is only the structural mass 
% Task:
%   -Estimate t to match results for nu = 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task1.possible_t = linspace(0,0.01,data.ntries);

task1.firstWB = 3.47;
task1.secondWB = 21.27;
task1.thirdWB = 58.13;
task1.chordwise = 21.75;
task1.torsion = 22.57;

% Cross sectional area of the beam as function of thickness [m^2]
task1.A = @(t) (data.nu*data.Cstr^2)-((data.nu*data.Cstr-2*data.tskin(t))*(data.Cstr-2*data.tspar(t)));
% Beam total mass as function of thickness [kg]
task1.Mass = @(t) data.rho*data.L*task1.A(t);
% Second moments of area of the beam as function of thickness [mm^4]
task1.Ixx = @(t) (((1/12)*data.Cstr*(data.nu*data.Cstr)^3)-((1/12)*(data.Cstr-2*data.tspar(t))*(data.nu*data.Cstr-2*data.tskin(t))^3));
task1.Izz = @(t) (((1/12)*data.nu*data.Cstr*(data.Cstr)^3)-((1/12)*(data.nu*data.Cstr-2*data.tskin(t))*(data.Cstr-2*data.tspar(t))^3));
task1.Iyy = @(t) ((1/12)*data.nu*data.Cstr*data.Cstr*(data.Cstr^2 + (data.nu*data.Cstr)^2)) - ((1/12)*(data.nu*data.Cstr-2*data.tskin(t))*(data.Cstr-2*data.tspar(t))*((data.Cstr-2*data.tspar(t))^2 + (data.nu*data.Cstr-2*data.tskin(t))^2));
task1.J   = @(t) data.rho*task1.Iyy(t);
% Vector with possible thickness [m]

% for i=1:length(task1.possible_t)
%     %Initializing the matrices
%     task1.K = zeros(data.nodes*6,data.nodes*6);
%     task1.masses = zeros(1,data.nodes);
%     task1.M = zeros(data.nodes*6,data.nodes*6);
%     task1.M_minus = zeros(data.nodes*6,data.nodes*6);
% 
%     current_index = 1;
% 
%     for j=1:(data.nodes-1)
%         task1.K(current_index:current_index+11,current_index:current_index+11) = ...
%             task1.K(current_index:current_index+11,current_index:current_index+11) + ...
%                     Stiffness_matrix_beam(data,task1.A(task1.possible_t(i)), ...
%                                           task1.Ixx(task1.possible_t(i)), ...
%                                           task1.Izz(task1.possible_t(i)), ...
%                                           task1.J(task1.possible_t(i)));
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
%         current_index = current_index + 6;
%     end
% 
%     current_index = 1;
% 
%     for j=1:data.nodes
%         task1.M(current_index:current_index+5,current_index:current_index+5) = ...
%                     Mass_matrix_beam(task1.masses(j), ...
%                                      task1.Ixx(task1.possible_t(i)), ...
%                                      task1.Iyy(task1.possible_t(i)), ...
%                                      task1.Izz(task1.possible_t(i)));
%         current_index = current_index + 6;
%     end
%     
%     for j=1:(6*data.nodes)
%         task1.M_minus(j,j) = task1.M(j,j)^(-1/2);
%     end
% 
%     task1.K_changed = task1.M_minus*task1.K*task1.M_minus;
% 
%     [task1.eigenvectors,task1.eigenvalues] = eig(task1.K_changed);
% 
%     task1.freqs = diag(real(sqrt(real(task1.eigenvalues))));
%     task1.modes = task1.M_minus*real(task1.eigenvectors);
%     
%     % Call function to solve task 01
%     %[task1.thickness,task1.area] = task1Fcn(data,task1);
%     if (task1.thickness == 0) || (task1.area <=0)
%         fprintf('No convergence in task 01 \n')
%     else
%         fprintf('Convergence in task 01, thickness = %f m and cross sectional area = %f m \n',task1.thickness,task1.area)
%     end
% 
% end

%% TASK 02

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumptions: 
%   -Wing is a continuous beam clamped at the wing root
%   -The only contribution to stiffness is from the wing box
%   -The mass is only the structural mass 
% Task:
%   -Estimate nu to match results for table 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% Vector with possible nu
task2.possible_nu = linspace(0.01,0.5,100000);

% Call function to solve task 02
%[task2.nu,task2.area] = task2Fcn(data,task2);
% if (task2.nu == 0) || (task2.area <=0)
%     fprintf('No convergence in task 02 \n')
% else
%     fprintf('Convergence in task 02, nu = %f m and cross sectional area = %f m \n',task2.nu,task2.area)
% end

%% TASK 03

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
task4.Kt = 10000000;
task4.EA = 0.35*task4.c;
task4.RS = 0.55*task4.c;
task4.FS = 0.15*task4.c;
task4.AF = 0.25*task4.c;

% ISA data
data.T0 = 288.15; 
data.alpha = -0.0065;
data.rho0 = 1.225;
data.g = 9.80665;
data.R = 287.1;
data.gamma = 1.4;
data.a = @(h) sqrt(data.gamma*data.R*(data.T0+data.alpha*h));
data.rho = @(h) data.rho0*(1+(data.alpha*h)/data.T0)^(-1-(data.g/(data.R*data.alpha)));
data.T = @(h) data.T0 + data.alpha*h;

task4.rhoSL = data.rho(0);
task4.aSL = data.a(0);
task4.qdiv = task4.Kt/(task4.S*task4.e*task4.Cl_alpha);
task4.vSL = sqrt(task4.qdiv/(0.5*task4.rhoSL));
task4.MSL = task4.vSL/task4.aSL;
task4.qdivi = task4.qdiv/sqrt(1-task4.MSL^2);

task4.KTAS = task4.vSL*1.94384;

fprintf('The incompressible divergence dynamic pressure at sea level qdiv = %f \n', task4.qdivi)
fprintf('The divergence speed at sea level Vdiv = %f knots \n', task4.KTAS)


%% TASK 05

task5.hRange = [0 10000 20000 30000]*0.3048;

for i = 1:length(task5.hRange)
    task5.rho(i) = data.rho(task5.hRange(i));
    task5.a(i) = data.a(task5.hRange(i));
    task5.A = (0.25*(task5.rho(i)^2)*(task5.a(i)^4))/(task4.qdiv^2);
    task5.coefvct = [-task5.A  0 task5.A  0  1 0 -1];     % Coefficient Vector
    task5.x = real(roots(task5.coefvct));
    task5.param = 0;
    for j = 6:-1:1
        task5.a = task5.x(j)-1;
    if ((task5.x(j) > task5.param) && (abs(task5.a) > 1e-03))
        task5.M(i) = task5.x(j);
        break
    end
    end
end
figure(1)
plot(task5.M, task5.hRange, 'b')
for i = 1:length(task5.hRange)
    yline(task5.hRange(i), 'r');
end
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

