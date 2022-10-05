function M = Mass_matrix_beam(m,Ix,Iy,Iz,rho,l)
%Returns the mass matrix for a beam

M = zeros(2,2);
% M(1,1) = m;
% M(2,2) = rho*Iy;

M(1,1) = 140;
M(2,1) = 70;
M(2,2) = M(1,1);
M(1,2) = M(2,1);

M = rho*l*M/420;

% M = zeros(6,6);
% M(1,1) = 0;
% M(2,2) = m;
% M(3,3) = 0;
% M(4,4) = 0;
% M(5,5) = 0;
% M(6,6) = Iy;

end

