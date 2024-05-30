function [f_x,f_y] = sbfem(d_b,l_b,c,omega,mu,n_x,n_y,epsilon,...
        epsilon_dot,X_att,X_att_dot,epsilon_constr,n_tay,...
        Vd_mat,Ld_vec,flag_tay,red,pb1,pb2)

% AUTHOR, AFFILIATION, DATE
% Simon Pfeil, OvGU Magdeburg, 30.05.24

% DESCRIPTION
% Solution of the Reynolds equation using the SBFEM (Scaled Boundary Finite 
% Element Method) and, if desired, with Taylor series approximations of 
% the eigenvalues and eigenvectors involved in this SBFEM solution; this
% function is used by the function 'run_up.m'

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
% - epsilon_constr = relative eccentricity where the Taylor series' which
%   will be used in this SBFEM solution were constructed (i.e., where the
%   eigenvalue and eigenvector derivatives were computed) [-]
% - n_tay = order of the Taylor series (maximum power in the Taylor 
%   series) [-]
% - Vd_mat = matrix containing the eigenvector derivatives for the Taylor
%   series [-]
% - Ld_vec = vector containing the eigenvalue derivatives for the Taylor
%   series [-]
% - flag_tay = determines whether to solve the eigenvalue problem
%   (flag_tay=0) or to use the Taylor series approximations (flag_tay=1)
%   [-]
% - red = number of considered modes for modal reduction of pressure
%   field [-]
% - pb1, pb2 = pressures at the two bearing boundaries (zero corresponds 
%   to atmospheric pressure) [Pa]

% OUTPUT VARIABLES
% - f_x = horizontal bearing force acting on the shell and, in the
%   opposite direction, on the shaft [N]
% - f_y = vertical bearing force acting on the shell and, in the
%   opposite direction, on the shaft [N]



% assembly


u = omega*(d_b/2);                                                         % circumferential surface velocity
p_ref = u*mu*(d_b/2)/(2*c^2);                                              % reference pressure due to nondimensionalization

L_X = 2*pi/n_x;                                                            % angular circumferential sector side length
gamma = l_b/d_b;                                                           % slenderness ratio

X_vec = linspace(0,2*pi-L_X,n_x)';                                         % vector of nodal circumferential coordinate values

H_vec = 1-epsilon*cos(X_vec);                                              % vector of nondimensionalized nodal gap functions
H3_vec = H_vec.^3;                                                         % vector of cubed nondimensionalized nodal gap functions
H_dot_vec = (d_b/u)*(-epsilon_dot*cos(X_vec)...                            % rate of change of the nondimensionalized nodal gap functions
    -epsilon*X_att_dot*sin(X_vec));

index_E_vec = vertcat(linspace(2,n_x,n_x-1)',1);                           % node number of eastern neighboring node for every node
index_W_vec = vertcat(n_x,linspace(1,n_x-1,n_x-1)');                       % node number of western neighboring node for every node

if flag_tay == 0                                                           % if no Taylor series will be used, i.e., if the eigenvalue problem will be solved
    
    index_0_vec = linspace(1,n_x,n_x)';                                    % vector containing all node numbers
    coef_E_vec = 1/(24*L_X)*(H3_vec+H3_vec(index_E_vec,1));                % vector containing the coefficients of E2 describing the interaction with the eastern neighbor
    coef_W_vec = 1/(24*L_X)*(H3_vec+H3_vec(index_W_vec,1));                % vector containing the coefficients of E2 describing the interaction with the western neighbor
    E2_mat = diag( coef_E_vec+coef_W_vec );                                % construction of E2-matrix: interaction of every node with itself
    E2_mat(n_x*(index_E_vec-1)+index_0_vec) = -coef_E_vec;                 % construction of E2-matrix: interaction with the eastern neighbor
    E2_mat(n_x*(index_W_vec-1)+index_0_vec) = -coef_W_vec;                 % construction of E2-matrix: interaction with western neighbor
    
end

diagE0_vec = (L_X/(12*gamma^2))*H3_vec;                                    % diagonal of E0, which is a diagonal matrix due to mass lumping

if flag_tay == 0                                                           % if the eigenvalue problem will be solved (no Taylor series) 
    trafo_vec = diagE0_vec.^(-0.5);                                        % vector representing the diagonal matrix that transforms the eigenvalue problem
end

R_vec = 0.5*(H_vec(index_E_vec,1)-H_vec(index_W_vec,1)) + L_X*H_dot_vec;   % right-hand side vector



% compute eigenvalues and -vectors by Taylor series if desired


if flag_tay == 1                                                           % if usage of a Taylor series is desired
    
    Delta_eps = epsilon-epsilon_constr;                                    % difference in relative eccentricity between the current shaft position and the one where the Taylor series coefficients were derived
    
    V_mat = zeros(n_x,red);                                                % initialization of the matrix that will later contain the eigenvectors
    L_vec = zeros(red,1);                                                  % initialization of the vector that will later contain the squared eigenvalues
    
    for j = 0:n_tay                                                        % loop through all powers involved in the Taylor series 
        fac = Delta_eps^j/factorial(j);                                    % factor to multiply coefficient by
        V_mat = V_mat + Vd_mat(:,(1+j*red):((j+1)*red))*fac;               % add contribution of Taylor series coefficient corresponding to current power to the matrix of eigenvectors
        L_vec = L_vec + Ld_vec(1,(1+j*red):((j+1)*red))'*fac;              % add contribution of Taylor series coefficient corresponding to current power to the vector of eigenvalues
    end
    
    V_mat = V_mat.*(diag(V_mat'*(diagE0_vec.*V_mat)).^(-0.5))';            % normalization
    
end



% compute eigenvalues and -vectors directly if desired


if flag_tay == 0                                                           % if no Taylor series is used
    
    T_mat = trafo_vec.*E2_mat.*trafo_vec';                                 % compute matrix of standard eigenvalue problem (transformation from generalized to standard)
    T_mat = (T_mat+T_mat')/2;                                              % ensure that rounding errors don't make the matrix unsymmetric

    [V_all_mat,L_all_mat,~] = eig(T_mat);                                  % solve eigenvalue problem

    [L_all_vec,order_vec] = sort(diag(L_all_mat));                         % sort eigenvalues
    V_all_mat = V_all_mat(:,order_vec);                                    % sort eigenvectors
    
    V_mat = V_all_mat(:,1:red);                                            % reduce set of eigenvectors ...
    L_vec = L_all_vec(1:red,1);                                            % ... and eigenvalues
    
end



% compute particular solution and integration constants


if flag_tay == 0                                                           % if no Taylor series has been used (which means that the eigenvalue problem has been solved)
    V_mat = trafo_vec.*V_mat;                                              % transform eigenvectors so that they satisfy the generalized (as opposed to standard) eigenvalue problem
end

% P_par_vec = -(V_mat*(vertcat(0,1./L_vec(2:red,1)).*(V_mat'*R_vec)));     % particular solution
% C_vec = -(V_mat'*(diagE0_vec.*P_par_vec));                               % integration constants

C_vec = vertcat(0,1./L_vec(2:red,1)).*((V_mat')*R_vec);                    % integration constants
P_par_vec = -(V_mat*C_vec);                                                % particular solution



% compute bearing forces


lambda_vec = vertcat(0,sqrt(L_vec(2:red,1)));                              % vector of eigenvalues (L_vec actually contains the squared eigenvalues; hence, the square root)

if pb1 == 0 && pb2 == 0                                                    % if the pressures prescribed at the bearing boundaries are nonzero
    
    % analytical integration in the axial direction
    
    p_bar_vec = p_ref * ( P_par_vec + V_mat*(vertcat(1,tanh(...            % averaging of the pressures in the axial direction based on analytical integration
        lambda_vec(2:red,1))./lambda_vec(2:red,1)).*C_vec) );
    p_bar_vec = ( p_bar_vec + abs(p_bar_vec) ) / 2;                        % Guembel condition
    
else
    
    % numerical integration in the axial direction
    
    xi_vec = linspace(0,1,n_y);                                            % vector containing the axial data points for evaluation of the pressure field (corresponding to the axial nodal positions in the FVM) for one half of the bearing, using a dimensionless coordinate xi equal to 0 at the axial bearing center and 1 at the bearing boundary
    
    P_mat_half = (V_mat*diag(C_vec./cosh(lambda_vec)))*...                 % nondimensionalized pressure field for one half of the bearing
        cosh(lambda_vec*xi_vec) + P_par_vec;
    
    p_mat_half = P_mat_half*p_ref;                                         % conversion from dimensionless to physical pressure (not yet considering periodic node, second half of the bearing, and Guembel condition)
    
    p_mat = horzcat(fliplr(p_mat_half(:,2:n_y)),p_mat_half);
    p_mat = p_mat + ( (pb1+pb2)/2 + ...                                    % consideration of nonzero boundary pressures: the solution computed for the inhomogeneous differential equation with homogeneous BCs is superposed with a linear function which satisfies the homogeneous differential equation and the inhomogeneous BCs
        ((pb2-pb1)/2)*horzcat(-fliplr(xi_vec),xi_vec(1,2:n_y)) );
    p_mat = (p_mat+abs(p_mat))/2;                                          % Guembel condition
    
    av_y_vec = vertcat(0.5,ones(2*n_y-3,1),0.5)/(2*n_y-2);                 % weigths for summation of the pressures in the axial direction in order to compute the axial average
    p_bar_vec = p_mat*av_y_vec;                                            % averaging of the pressures in the axial direction
    
end

% numerical integration in the circumferential direction

f_x = (cos((X_vec+X_att)')*p_bar_vec)*(l_b*L_X*(d_b/2));               	   % summation/integration in the circumferential direction, leading to the hydrodynamic force acting on the shell (and, with opposite sign, on the shaft) in the horizontal direction
f_y = (sin((X_vec+X_att)')*p_bar_vec)*(l_b*L_X*(d_b/2));               	   % summation/integration in the circumferential direction, leading to the hydrodynamic force acting on the shell (and, with opposite sign, on the shaft) in the vertical direction



end
