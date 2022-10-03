%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOMEWORK 1 AEROELASTICITY
% STUDY OF DIVERGENCE
% ISAAC ROBLEDO MARTIN
% FERNANDO RUIZ CERRAJERO
% JUAN ALFARO MORENO
% UC3M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%% TASK 01

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumptions: 
%   -Wing is a continuous beam clamped at the wing root
%   -The only contribution to stiffness is from the wing box
%   -The mass is only the structural mass 
% Task:
%   -Estimate t to match results for nu = 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task1.firstWB = 3.47;
task1.secondWB = 21.27;
task1.thirdWB = 58.13;
task1.chordwise = 21.75;
task1.torsion = 22.57;

% Cross sectional area of the beam as function of thickness [m^2]
task1.A = @(t) (data.nu*data.Cstr^2)-((data.nu*data.Cstr-2*data.tskin(t))*(data.Cstr-2*data.tspar(t)));
% Second moments of area of the beam as function of thickness [mm^4]
task1.Ixx = @(t) (((1/12)*data.nu*data.Cstr^4)-((1/12)*(data.Cstr-2*data.tspar(t))*(data.nu*data.Cstr-2*data.tskin(t))^3));
task1.Izz = @(t) (((1/12)*data.nu*data.Cstr^4)-((1/12)*((data.Cstr-2*data.tspar(t))^3)*(data.nu*data.Cstr-2*data.tskin(t))));
% Vector with possible thickness [m]
task1.possible_t = linspace(0.01,1,100000);

% Call function to solve task 01
[task1.thickness,task1.area] = task1Fcn(data,task1);
if (task1.thickness == 0) || (task1.area <=0)
    fprintf('No convergence in task 01 \n')
else
    fprintf('Convergence in task 01, thickness = %f m and cross sectional area = %f m \n',task1.thickness,task1.area)
end


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
[task2.nu,task2.area] = task2Fcn(data,task2);
if (task2.nu == 0) || (task2.area <=0)
    fprintf('No convergence in task 02 \n')
else
    fprintf('Convergence in task 02, nu = %f m and cross sectional area = %f m \n',task2.nu,task2.area)
end

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

task4.S = 0.625; 
task4.S = 62.5;
task4.c = 6.25; 
task4.Cl_alpha = 2*pi;
task4.Kt = 1;
task4.EA = 0.35*task4.c;
task4.RS = 0.55*task4.c;
task4.FS = 0.15*task4.c;
task4.AF = 0.25*task4.c;

% ISA data
task4.T0 = 288.15; 
task4.alpha = -0.0065;
task4.rho0 = 1.225;
task4.g = 9.80665;
task4.R = 287.1;

% Dynamic pressure computation
task4.qdiv = task4.Kt/(task4.Cl_alpha*task4.S*task4.e);

task4.hVector = linspace(0,11000,1000);

for i = 1:length(task4.hVector)
    task4.vVector(i) = sqrt(1/(2*task4.qdiv*task4.rho0*(1+(task4.alpha*task4.hVector(i))/task4.T0)^(-1 ...
        -task4.g/(task4.R*task4.alpha))));
end

plot(task4.vVector,task4.hVector)

%% FUNCTIONS
function [t,area] = task1Fcn(data,task1)
    for i = 1:length(task1.possible_t)
         area = task1.A(task1.possible_t(i));
         C1 = (data.rho*area*data.L*task1.firstWB^2)/(data.E*task1.Izz(task1.possible_t(i)));
         C2 = (data.rho*area*data.L*task1.secondWB^2)/(data.E*task1.Izz(task1.possible_t(i)));
         C3 = (data.rho*area*data.L*task1.secondWB^2)/(data.E*task1.Izz(task1.possible_t(i)));
         
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

