% AUTHOR, AFFILIATION, DATE
% Simon Pfeil, OvGU Magdeburg, 30.05.24

% DESCRIPTION
% Computation of the eigenvalue/eigenvector derivatives of an SBFEM
% simulation model for plain bearings; these derivatives can then be
% used as Taylor series coefficients for approximating the eigenvalues
% and eigenvectors of the SBFEM solution during a time integration, which
% avoids solving the eigenvalue problem at every time step; this script
% uses the function 'construct_series.m'; run this script in Matlab



% clear variables, clear console, close figures, etc. ...


clear variables
close all
clc
addpath(genpath(pwd()))
dbstop if error



% parameters


n_x = 100;                                                                 % circumferential number of nodes
n_tay = 2;                                                                 % maximum power in the Taylor series'
tol_eignvl = 1e-4;                                                         % if the relative difference between two eigenvalues is less than this tolerance, they are considered identical
gamma_ref = 1;                                                             % arbitrarily chosen slenderness ratio (bearing length over bearing diameter); the generated data can be adjusted to any other slenderness ratio later



% compute and save eigenvalue/eigenvector derivatives


tic0 = tic;                                                                % start clock for computational time

load('1_locations.mat','constr_vec')                                       % constr_vec contains the points throughout the range of relative eccentricities where Taylor series' will be constructed and ...
load('2_reductions.mat','red_vec')                                         % ... red_vec describes, for each of these points, the size of the subset of modes to be considered in the computation of the pressure field ("red" for reduction, as in modal reduction)

n_constr = length(constr_vec);                                             % number of points where Taylor series' will be constructed
n_col = sum(red_vec*(n_tay+1));                                            % number of columns of the matrices for storing the eigenvalue and eigenvector derivatives

Vd_allpoints_mat = zeros(n_x,n_col);                                       % matrix for storing all eigenvectors and their derivatives at all data points
Ld_allpoints_vec = zeros(1,n_col);                                         % matrix for storing all eigenvalues and their derivatives at all data points

ind = 1;                                                                   % variable for trackig the first column index corresponding to the eigenvalue/eigenvector derivatives at a given data point 
ind_vec = zeros(n_constr,1);                                               % vector for storing the first column index corresponding to the eigenvalue/eigenvector derivatives at a given data point, for every data point 
for i = 1:n_constr                                                         % loop through all points
    epsilon_constr = constr_vec(i,1);                                      % current point
    red = red_vec(i,1);                                                    % number of modes to consider at the current point
    [Vd,Ld] = construct_series(gamma_ref,epsilon_constr,...                % compute eigenvectors and eigenvalues and their derivatives at current point
        n_x,n_tay,tol_eignvl);
    for j = 1:(n_tay+1)                                                    % loop through all derivative orders
        Vd_allpoints_mat(:,(ind+(j-1)*red):(ind-1+j*red)) = Vd(:,1:red,j); % store j-th derivative of the eigenvectors 
        Ld_allpoints_vec(1,(ind+(j-1)*red):(ind-1+j*red)) = Ld(1:red,j)';  % store j-th derivative of the eigenvalues
    end
    ind_vec(i,1) = ind;                                                    % store column number where the data corresponding to the current point begins
    ind = ind + red*(n_tay+1);                                             % calculate the column number where the data corresponding to the next point will begin
end

save('3_coefficients.mat','n_x','n_tay','gamma_ref',...                    % save results and parameters
    'Vd_allpoints_mat','Ld_allpoints_vec','ind_vec')

t_pre = toc(tic0);                                                         % stop clock for compuational time

disp(['Elapsed time is ',num2str(t_pre),' seconds.'])

