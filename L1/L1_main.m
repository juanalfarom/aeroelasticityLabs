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
    print('No convergence in task 01')
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

