function [f_x,f_y] = fvm(d_b,l_b,c,omega,mu,n_x,n_y,epsilon,...
        epsilon_dot,X_att,X_att_dot,pb1,pb2)

% AUTHOR, AFFILIATION, DATE
% Simon Pfeil, OvGU Magdeburg, 30.05.24

% DESCRIPTION
% Solution of the Reynolds equation using the FVM (Finite Volume Method);
% this function is used by the function 'run_up.m'

% INPUT VARIABLES
% - d_b = bearing diameter [m]
% - l_b = bearing length [m]
% - c = radial clearance [m]
% - omega = angular velocity of the shaft [rad/s]
% - mu = dynamic viscosity [Pa*s]
% - n_x = circumferential number of nodes [-]
% - n_y = axial number of nodes [-]
% - epsilon = relative eccentricity [-]
% - epsilon_dot = rate of change of the relative eccentricity [1/s]
% - X_att = attitude angle [rad]
% - X_att_dot = rate of change of attitude angle [rad/s]
% - pb1, pb2 = pressures at the two bearing boundaries (zero corresponds 
%   to atmospheric pressure) [Pa]

% OUTPUT VARIABLES
% - f_x = horizontal bearing force acting on the shell and, in the
%   opposite direction, on the shaft [N]
% - f_y = vertical bearing force acting on the shell and, in the
%   opposite direction, on the shaft [N]



% some preperations


dof_y = n_y-1;                                                             % number of DOFs in the axial direction (axial number of nodes minus 1 because the nodes at the bearing boundary have Dirichlet BCs and therefore aren't DOFs)
dof = n_x*dof_y;                                                           % total number of DOFs
L_X = 2*pi/n_x;                                                            % dimensionless circumferential control volume side length
L_Y = (l_b/d_b)/(n_y-1);                                                   % dimensionless axial control volume side length
p_ref = 6*(d_b/2)^2*mu*omega/c^2;                                          % reference pressure according to nondimensionalization
X_vec = linspace(0,2*pi-L_X,n_x)';                                         % circumferential nodal coordinates
H_vec = 1-epsilon*cos(X_vec);                                              % dimensionless nodal gap widths (which only vary in the circumferential direction, since shaft tilting is neglected)
H3_vec = H_vec.^3;                                                         % cubed dimensionless gap widths
H_dot_vec = -epsilon_dot*cos(X_vec)-epsilon*X_att_dot*sin(X_vec);          % rate of change of the dimensionless gap widths



% determine DOF numbers and respective neighboring DOF numbers


center = (1:n_x)';                                                         % DOFs at the bearing center
other = (n_x+1:n_x*dof_y)';                                                % DOFs not at the bearing center
index_0_vec = (1:dof)';                                                    % all DOFs
index_N_vec = vertcat(center+n_x,other+n_x);                               % neighboring DOF number in the northern direction for every DOF (the last n_x entries are at the bearing boundary, but no DOFs exist there, so these entries will not be referenced)
index_S_vec = vertcat(center+n_x,other-n_x);                               % neighboring DOF number in the southern direction for every DOF (the first n_x entries are replaced by the northern DOFs, enforcing the symmetric BC at the bearing center)
index_E_vec = index_0_vec + 1 - sparse(n_x*(1:dof_y)',1,n_x,dof,1);        % neighboring DOF number in the eastern direction for every DOF
index_W_vec = index_0_vec - 1 + sparse(1+n_x*(0:dof_y-1)',1,n_x,dof,1);    % neighboring DOF number in the western direction for every DOF



% compute coefficients on the LHS of the discretized Reynolds equation


a_SN_vec = repmat((L_X/L_Y)*H3_vec,dof_y,1);                               % the southern and northern coefficients are identical because the gap function is assumed to be constant in the axial direction
a_E_vec = repmat((L_Y/(2*L_X))*...                                         % eastern coefficients
    (H3_vec+H3_vec(index_E_vec(1:n_x,1),1)),dof_y,1);
a_W_vec = repmat(a_E_vec(index_W_vec(1:n_x,1),1),dof_y,1);                 % western coefficients



% construct matrix under consideration of the BCs


rows_vec = vertcat(index_0_vec,index_0_vec,index_0_vec,index_0_vec,...     % row numbers of the nonzero matrix entries stored in a vector
    index_0_vec(1:(dof-n_x),1));
cols_vec = vertcat(index_0_vec,index_W_vec,index_E_vec,index_S_vec,...     % column numbers of the nonzero matrix entries stored in a vector
    index_N_vec(1:(dof-n_x),1));
coef_vec = vertcat(-a_W_vec-a_E_vec-2*a_SN_vec,a_W_vec,a_E_vec,...         % nonzero matrix entries stored in a vector
    a_SN_vec,a_SN_vec(1:(dof-n_x),1));
entries_center_vec = vertcat(center,center+dof,center+2*dof,...            % this vector points to those entries in coef_vec that evaluate the Reynolds equation at the bearing center; these entries ...
    center+3*dof,center+4*dof);
coef_vec(entries_center_vec,1) = coef_vec(entries_center_vec,1)/2;         % ... must be halved because of the half-sized control volumes (otherwise the matrix would be unsymmetric)
K_mat = sparse(rows_vec,cols_vec,coef_vec,dof,dof);                        % construction of the matrix



% construct RHS vector


R_vec = repmat((L_Y/2)*(H_vec(index_E_vec(1:n_x,1),1)-H_vec(...            % construction of the RHS vector
    index_W_vec(1:n_x,1),1))+((2/omega)*L_X*L_Y)*H_dot_vec,dof_y,1);
R_vec(1:n_x,1) = R_vec(1:n_x,1)/2;                                         % the entries at the bearing center are halved because of the half-sized control volumes



% computation of the dimensionless pressures at all DOFs


P_vec_dof = K_mat\R_vec;                                                   % the linear system of equations is solved; the nodal solutions are stored in a vector



% postprocessing of the computed pressure distribution


P_vec = vertcat(P_vec_dof,zeros(n_x,1));                                   % inclusion of the pressures at the bearing boundary, which are equal to zero (the modifications necessary for the case of inhomogeneous BCs are performed below)
P_mat_half = reshape(P_vec,[n_x,n_y]);                                     % the pressure field in the considered half of the bearing is stored in a matrix
P_mat = horzcat(fliplr(P_mat_half(:,2:n_y)),P_mat_half);                   % the whole pressure field is stored in a matrix (the considered half is mirrored)
p_mat = P_mat*p_ref;                                                       % the dimensionless pressures are converted to physical pressures in [Pa]
if pb1 ~= 0 || pb2 ~= 0                                                    % if the boundary conditions are inhomogeneous
    xi_vec = linspace(0,1,n_y);                                            % a vector containing the axial nodal positions normalized in such a way that -1 and 1 correspond to the bearing boundaries
    p_mat = p_mat + ( (pb1+pb2)/2 + ...                                    % the solution computed for the inhomogeneous differential equation with homogeneous BCs is superposed with a linear function which satisfies the homogeneous differential equation and the inhomogeneous BCs
        ((pb2-pb1)/2)*horzcat(-fliplr(xi_vec),xi_vec(1,2:n_y)) );
end
p_mat = (p_mat+abs(p_mat))/2;                                              % Guembel condition



% computation of the bearing forces by integration of the pressure


pdy = p_mat*vertcat(0.5,ones(2*n_y-3,1),0.5);                              % summation of the nodal pressures in the axial direction (the ones at the boundary are weighted only with 0.5)
f_x = ( ones(1,n_x)*(pdy.*cos(X_vec+X_att)) )*(L_X*L_Y*(d_b^2/4));         % horizontal component of hydrodynamic force acting on the shell (and, with opposite sign, on the shaft)
f_y = ( ones(1,n_x)*(pdy.*sin(X_vec+X_att)) )*(L_X*L_Y*(d_b^2/4));         % vertical component of hydrodynamic force acting on the shell (and, with opposite sign, on the shaft)



end