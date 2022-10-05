function K = Stiffness_matrix_beam(data,A,Ix,Iz,J,L)
%Calculate the stiffness matrix for a given beam

K = zeros(4,4);

K(1,:) = [12*data.E*Ix/(L^3),0,-12*data.E*Ix/(L^3),0];
K(3,:) = -K(1,:);
K(2,:) = [0,data.G*J/L,0,-data.G*J/L];
K(4,:) = -K(2,:);


% K = zeros(12,12);
% 
% K(1,:) = [12*data.E*Iz/(data.L^3),0,0,0,0,-6*data.E*Iz/(data.L^2),-12*data.E*Iz/(data.L^3),0,0,0,0,-6*data.E*Iz/(data.L^2)];
% K(2,:) = [0,data.E*A/data.L,0,0,0,0,0,-data.E*A/data.L,0,0,0,0];
% K(3,:) = [0,0,12*data.E*Ix/(data.L^3),-6*data.E*Ix/(data.L^2),0,0,0,0,-12*data.E*Ix/(data.L^3),-6*data.E*Ix/(data.L^2),0,0];
% K(4,:) = [0,0,-6*data.E*Ix/(data.L^2),4*data.E*Ix/(data.L),0,0,0,0,6*data.E*Ix/(data.L^2),2*data.E*Ix/(data.L),0,0];
% K(5,:) = [0,0,0,0,data.G*J/data.L,0,0,0,0,0,-data.G*J/data.L,0];
% K(6,:) = [-6*data.E*Iz/(data.L^2),0,0,0,0,4*data.E*Iz/(data.L),6*data.E*Iz/(data.L^2),0,0,0,0,2*data.E*Iz/(data.L)];
% K(7,:) = [-12*data.E*Iz/(data.L^3),0,0,0,0,6*data.E*Iz/(data.L^2),12*data.E*Iz/(data.L^3),0,0,0,0,6*data.E*Iz/(data.L^2)];
% K(8,:)= -K(2,:);
% K(9,:)= -K(3,:);
% K(10,:) = [0,0,-6*data.E*Ix/(data.L^2),2*data.E*Ix/(data.L),0,0,0,0,6*data.E*Ix/(data.L^2),4*data.E*Ix/(data.L),0,0];
% K(11,:)= -K(5,:);
% K(12,:) = [-6*data.E*Iz/(data.L^2),0,0,0,0,2*data.E*Iz/(data.L),6*data.E*Iz/(data.L^2),0,0,0,0,4*data.E*Iz/(data.L)];

end

