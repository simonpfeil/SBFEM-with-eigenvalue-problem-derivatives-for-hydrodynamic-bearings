function [Vd_mats,Ld_mat] = construct_series(gamma,epsilon_constr,...
    n_x,n_tay,tol_eignvl)

% AUTHOR, AFFILIATION, DATE
% Simon Pfeil, OvGU Magdeburg, 30.05.24

% DESCRIPTION
% Computation of the Taylor series coefficients for approximation of the 
% eigenvalues and eigenvectors involved in the SBFEM solution of the 
% Reynolds equation; the independent variable of these Taylor series' is 
% the relative eccentricity epsilon; this function is used by the Matlab
% script 'A_precomputation.m'

% INPUT VARIABLES
% - gamma = slenderness ratio [-]
% - epsilon_constr = the relative eccentricity where the Taylor series 
%   will be constructed (i.e., where the approximation will be exact) [-]
% - n_x = circumferential number of nodes [-]
% - n_tay = order of the Taylor series (for example, n_tay=2 means that 
%   the highest power in the Taylor series is 2) [-]
% - tol_eignvl = if the difference between two eigenvalues is below this 
%   prescribed tolerance, these two eigenvalues will be treated as 
%   identical [-]

% OUTPUT VARIABLES
% - Vd_mats = array containing the computed eigenvector derivatives [-]
% - Ld_mat = array containing the computed eigenvalue derivatives [-]



% assembly


L_X = 2*pi/n_x;                                                            % angular circumferential sector side length

X_vec = linspace(0,2*pi-L_X,n_x)';                                         % angular circumferential positions of the nodes
cos_vec = cos(X_vec);                                                      % cosines of the angular circumferential positions of the nodes
H_vec = 1-epsilon_constr*cos_vec;                                          % nondimensionalized gap function at the nodes

H3_vec = H_vec.^3;                                                         % cubed nondimensionalized gap function at the nodes
Hd1_vec = -cos_vec;                                                        % derivative of the nodal nondimensionalized gap functions with respect to the relative eccentricity epsilon
H3d1_vec = 3*(H_vec.^2).*Hd1_vec;                                          % derivative of the cubed nodal nondimensionalized gap functions with respect to the relative eccentricity epsilon
H3d2_vec = 6*H_vec.*Hd1_vec.^2;                                            % 2nd derivative of the cubed nodal nondimensionalized gap functions with respect to the relative eccentricity epsilon
H3d3_vec = 6*Hd1_vec.^3;                                                   % 3rd derivative of the cubed nodal nondimensionalized gap functions with respect to the relative eccentricity epsilon

index_0_vec = linspace(1,n_x,n_x)';                                        % vector containing all node numbers
index_E_vec = vertcat(linspace(2,n_x,n_x-1)',1);                           % vector containing the node numbers of the eastern neighbors
index_W_vec = vertcat(n_x,linspace(1,n_x-1,n_x-1)');                       % vector containing the node numbers of the western neighbors

E2d_mats = zeros(n_x,n_x,min(3,n_tay));                                    % array for storing 1st, 2nd, and 3rd derivatives of the matrix E2 with respect to the eccentricity epsilon
E0d_mat = zeros(n_x,min(3,n_tay));                                         % array for storing 1st, 2nd, and 3rd derivatives of the matrix E0 (given by vectors representing the main diagonals) with respect to the eccentricity epsilon

if n_tay >= 1
    coef_E_vec = 1/(24*L_X)*(H3d1_vec+H3d1_vec(index_E_vec,1));            % coefficients of E2d1 (derivative of E2 with respect to the eccentricity epsilon) describing the nodal interaction with the eastern neighbor
    coef_W_vec = 1/(24*L_X)*(H3d1_vec+H3d1_vec(index_W_vec,1));            % coefficients of E2d1 (derivative of E2 with respect to the eccentricity epsilon) describing the nodal interaction with the western neighbor
    E2_mat = diag( coef_E_vec+coef_W_vec );                                % construction of the matrix E2d1 (derivative of E2 with respect to the eccentricity epsilon): main diagonal
    E2_mat(n_x*(index_E_vec-1)+index_0_vec) = -coef_E_vec;                 % construction of the matrix E2d1 (derivative of E2 with respect to the eccentricity epsilon): entries representing the nodal interaction with the eastern neighbor
    E2_mat(n_x*(index_W_vec-1)+index_0_vec) = -coef_W_vec;                 % construction of the matrix E2d1 (derivative of E2 with respect to the eccentricity epsilon): entries representing the nodal interaction with the western neighbor
    E2d_mats(:,:,1) = E2_mat;                                              % store E2d1
    E0d_mat(:,1) = (L_X/(12*gamma^2))*H3d1_vec;                            % construction of the matrix E0d1 (derivative of E0 with respect to the eccentricity epsilon), but in the form of a vector representing the main diagonal (E0d1 is a diagonal matrix)
end

if n_tay >= 2
    coef_E_vec = 1/(24*L_X)*(H3d2_vec+H3d2_vec(index_E_vec,1));            % coefficients of E2d2 (second derivative of E2 with respect to the eccentricity epsilon) describing the nodal interaction with the eastern neighbor
    coef_W_vec = 1/(24*L_X)*(H3d2_vec+H3d2_vec(index_W_vec,1));            % coefficients of E2d2 (second derivative of E2 with respect to the eccentricity epsilon) describing the nodal interaction with the western neighbor
    E2_mat = diag( coef_E_vec+coef_W_vec );                                % construction of the matrix E2d2 (second derivative of E2 with respect to the eccentricity epsilon): main diagonal
    E2_mat(n_x*(index_E_vec-1)+index_0_vec) = -coef_E_vec;                 % construction of the matrix E2d2 (second derivative of E2 with respect to the eccentricity epsilon): entries representing the nodal interaction with the eastern neighbor
    E2_mat(n_x*(index_W_vec-1)+index_0_vec) = -coef_W_vec;                 % construction of the matrix E2d2 (second derivative of E2 with respect to the eccentricity epsilon): entries representing the nodal interaction with the western neighbor
    E2d_mats(:,:,2) = E2_mat;                                              % store E2d2
    E0d_mat(:,2) = (L_X/(12*gamma^2))*H3d2_vec;                            % construction of the matrix E0d2 (second derivative of E0 with respect to the eccentricity epsilon), but in the form of a vector representing the main diagonal (E0d2 is a diagonal matrix)
end

if n_tay >= 3
    coef_E_vec = 1/(24*L_X)*(H3d3_vec+H3d3_vec(index_E_vec,1));            % coefficients of E2d3 (second derivative of E2 with respect to the eccentricity epsilon) describing the nodal interaction with the eastern neighbor
    coef_W_vec = 1/(24*L_X)*(H3d3_vec+H3d3_vec(index_W_vec,1));            % coefficients of E2d3 (second derivative of E2 with respect to the eccentricity epsilon) describing the nodal interaction with the western neighbor
    E2_mat = diag( coef_E_vec+coef_W_vec );                                % construction of the matrix E2d3 (third derivative of E2 with respect to the eccentricity epsilon): main diagonal
    E2_mat(n_x*(index_E_vec-1)+index_0_vec) = -coef_E_vec;                 % construction of the matrix E2d3 (third derivative of E2 with respect to the eccentricity epsilon): entries representing the nodal interaction with the eastern neighbor
    E2_mat(n_x*(index_W_vec-1)+index_0_vec) = -coef_W_vec;                 % construction of the matrix E2d3 (third derivative of E2 with respect to the eccentricity epsilon): entries representing the nodal interaction with the western neighbor
    E2d_mats(:,:,3) = E2_mat;                                              % store E2d3
    E0d_mat(:,3) = (L_X/(12*gamma^2))*H3d3_vec;                            % construction of the matrix E0d3 (third derivative of E0 with respect to the eccentricity epsilon), but in the form of a vector representing the main diagonal (E0d3 is a diagonal matrix)
end

coef_E_vec = 1/(24*L_X)*(H3_vec+H3_vec(index_E_vec,1));                    % coefficients of E2 describing the nodal interaction with the eastern neighbor
coef_W_vec = 1/(24*L_X)*(H3_vec+H3_vec(index_W_vec,1));                    % coefficients of E2 describing the nodal interaction with the western neighbor
E2_mat = diag( coef_E_vec+coef_W_vec );                                    % construction of the matrix E2: main diagonal
E2_mat(n_x*(index_E_vec-1)+index_0_vec) = -coef_E_vec;                     % construction of the matrix E2: entries representing the nodal interaction with the eastern neighbor
E2_mat(n_x*(index_W_vec-1)+index_0_vec) = -coef_W_vec;                     % construction of the matrix E2: entries representing the nodal interaction with the western neighbor

trafo_vec = ((L_X/(12*gamma^2))*H3_vec).^(-0.5);                           % this vector represents the diagonal matrix for transforming the eigenvalue from generalized to standard
T_mat = trafo_vec.*E2_mat.*trafo_vec';                                     % matrix of standard eigenvalue problem



% solution of the eigenvalue problem


T_mat = (T_mat+T_mat')/2;                                                  % this operation should recover the symmetry of this matrix in case rounding errors have made it unsymmetric (not sure this actually happens, probably not)
[V_mat,L_mat,~] = eig(T_mat);                                              % solution of eigenvalue problem, yielding the 0th-order coefficients of the Taylor series; note that the generalized eigenvalue problem was transformed to standard (the generalized formulation will only be used for computation of the orders >0)

[L_vec,order_vec] = sort(diag(L_mat));                                     % the eigenvalues, given by the main diagonal entries of L_mat, are sorted in ascending order and stored in a vector L_vec; the changes in order are tracked by order_vec
V_mat = V_mat(:,order_vec);                                                % the eigenvectors, given by the columns of V_mat, are sorted so that the order matches the eigenvalues

V_mat = trafo_vec.*V_mat;                                                  % transformation of the eigenvectors so that the generalized eigenvalue problem (instead of the standard eigenvalue problem) is satisfied



% computation of the eigenvalue and eigenvector derivatives


i_vec = zeros(1,n_x);                                                      % vector for storing the respective index of the maximum absolute value within each eigenvector
fac_vec = zeros(1,n_x);                                                    % vector for storing some factors which will be used later during the definition of the main diagonal entries of Q

for j = 1:n_x                                                              % loop through all eigenvectors
    [~,i] = max(abs(V_mat(:,j)));                                          % determine the index i pointing to the maximum absolute value of the current eigenvector ...
    i_vec(1,j) = i;                                                        % ... and store this index for later use
    fac_vec(1,j) = -1/V_mat(i,j);                                          % the negative reciprocal of the absolute-wise largest entry of the eigenvector is stored for later use (for computation of the main diagonal entry of Q)
end

divide_delta_L_mat = zeros(n_x,n_x);                                       % matrix for storing the reciprocal of the difference between two eigenvalues, for all combinations of eigenvalues

for i = 2:n_x                                                              % loop through eigenvalues
    for j = 1:(i-1)                                                        % for each eigenvalue, loop through the other eigenvalues
        if abs((L_vec(i,1)-L_vec(j,1)))/((L_vec(i,1)+L_vec(j,1))/2) ...    % if the relative difference between the two eigenvalues is larger than the tolerance for considering the eigenvalues identical
                > tol_eignvl
            divide_delta_L = 1 / (L_vec(j,1)-L_vec(i,1));                  % compute reciprocal of difference between eigenvalues, ...
            divide_delta_L_mat(i,j) = divide_delta_L;                      % ... store this value,
            divide_delta_L_mat(j,i) = -divide_delta_L;                     % ... and store the negative counterpart of this value (for the case where the two eigenvalues are swapped)
        end
    end
end

Ld_mat = zeros(n_x,n_tay+1);                                               % array for storing all eigenvalue derivatives
Q_mats = zeros(n_x,n_x,n_tay+1);                                           % array for storing the Q matrices corresponding to the eigenvector derivatives
Vd_mats = zeros(n_x,n_x,n_tay+1);                                          % array for storing the eigenvector derivatives

Ld_mat(:,1) = L_vec;                                                       % store eigenvalues as 0th-order derivatives of the eigenvalues
Q_mats(:,:,1) = eye(n_x);                                                  % the Q matrix corresponding to the eigenvectors (i.e., to the 0th-order derivatives of the eigenvectors) is the identity matrix
Vd_mats(:,:,1) = V_mat;                                                    % store eigenvectors as 0th-order derivatives of the eigenvectors

for n = 1:n_tay                                                            % loop through all orders of derivative that need to be computed, from 1 up to the prescribed order of the Taylor series

    Kn = zeros(n_x,n_x);                                                   % the sum of all known terms of the n-th derivative of the eigenvalue problem will be stored in this matrix

    for k = 1:min(n,3)                                                     % loop through orders of derivative of the E2 and E0 matrices
        Kk = zeros(n_x,n_x);                                               % the sum of all terms related to the k-the derivative of E0 will be stored in this matrix
        for l = 0:(n-k)                                                    % loop through eigenvector derivatives (those of lower order than n, which are already known)
            Kk = Kk + Vd_mats(:,:,l+1).*(Ld_mat(:,n-k-l+1)'*...            % consider term with l-th eigenvector derivative
                (factorial(n)/(factorial(k)*factorial(l)*...
                factorial(n-k-l))));
        end
        Kn = Kn - (E0d_mat(:,k)).*Kk + ...                                 % consider terms related to the k-th derivative of E2 and E0
            nchoosek(n,k)*E2d_mats(:,:,k)*Vd_mats(:,:,n-k+1);
    end

    Kn = V_mat' * Kn;                                                      % some tranformation simplifying the equation

    for l = 1:(n-1)                                                        % loop through known Q matrices
        Kn = Kn - Q_mats(:,:,l+1).*(Ld_mat(:,n-l+1)'*nchoosek(n,l));       % consider term involving l-th Q matrix and not explicitly depending on E0 or E2
    end

    Ld_mat(:,n+1) = diag(Kn);                                              % store n-th eigenvalue derivative (given by the diagonal entries of Kn)
    Q_mats(:,:,n+1) = Kn.*divide_delta_L_mat;                              % compute and store n-th Q matrix

    for j = 1:n_x                                                          % loop through main diagonal elements of the Q matrix
        Q_mats(j,j,n+1) = fac_vec(1,j) * ...                               % employ convention for choosing the main diagonal element, which would not be uniquely determined otherwise (this has something to do with eigenvectors being scalable)
            (V_mat(i_vec(1,j),:)*Q_mats(:,j,n+1));
    end

    Vd_mats(:,:,n+1) = V_mat*Q_mats(:,:,n+1);                              % compute and store n-th eigenvector derivative

end



end
