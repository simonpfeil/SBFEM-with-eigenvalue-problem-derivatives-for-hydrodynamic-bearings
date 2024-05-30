function [t,z,n_calls] = run_up(d_b,l_b,c,mu,m,omegaspan,t_start,t_end,...
    z_0,method,g,unb,fac_frq,phi_u,pb1,pb2,mbf,n_x_noTaylor)

% AUTHOR, AFFILIATION, DATE
% Simon Pfeil, OvGU Magdeburg, 30.05.24

% DESCRIPTION
% Run-up simulation with a rotor represented by a point mass in a plain 
% bearing; the Reynolds equation is solved with either the FVM or the 
% SBFEM; this function is called by 'main.m' and it calls the ODE solver
% for the time integration, as well as the functions 'fvm.m' or 'sbfem.m' 
% to solve the Reynolds equation with the FVM or SBFEM, respectively; it 
% also contains the nested function 'odeFun' describing the equation of 
% motion

% INPUT VARIABLES
% - d_b = bearing diameter [m]
% - l_b = bearing length [m]
% - c = radial clearance [m]
% - mu = dynamic viscosity [Pa*s]
% - m = rotor mass [kg]
% - omegaspan = row vector with two entries describing the angular 
%   velocitiy of the shaft at the beginning and at the end, respectively, 
%   in [rad/s]; note that this function assumes a constant angular accele-
%   ration d(omega)/d(t)=(omegaspan(1,2)-omegaspan(1,1))/(t_end-t_start)
% - t_start = time at the beginning of the run-up [s]
% - t_end = time at the end of the run-up [s]
% - z_0 = column vector containing 4 entries representing the initial
%   conditions at t_start, namely the horizontal displacement [m], 
%   vertical displacement [m], horizontal velocity [m/s], and vertical 
%   velocity [m/s] of the shaft (in this order)
% - method = determines the method for solving the Reynolds equation:
%   1 - FVM, 2 - SBFEM with solution of eigenvalue problem at every call of
%   the Reynolds equation, 3 - SBFEM with Taylor series approximations of 
%   the eigenvalues and eigenvectors [-]
% - g = gravitation constant [N/kg]
% - unb = shaft unbalance [kg*m]
% - fac_frq = output frequency for the simulation results relative to the
%   maximum frequency of the shaft rotation (for example, fac_frq=5 in
%   combination with a fastest rotation of 1000Hz will lead to an output
%   frequency of 5000Hz, i.e., 5000 data points in the time domain for
%   every second) [-]
% - phi_u = angular position of the unbalance at t=t_start, where phi_u=0
%   corresponds to the direction in which the horizontal shaft displacement
%   is measured and phi_u=pi/2 corresponds to the direction in which the
%   vertical shaft displacement is measured (phi_u=pi and phi_u=3*pi/2
%   correspond to the negative counterparts of these directions) [rad]
% - pb1, pb2 = pressures at the two bearing boundaries (zero corresponds 
%   to atmospheric pressure) [Pa]
% - mbf = multiplyer for bearing forces [-]; this program assumes a system
%   with only a single bearing (it's more of an academic example than a
%   realistic technical system), but we can mimic the existence of
%   multiple identical bearings by multiplying the bearing forces by a
%   factor mbf equal to the number of bearings (this may be useful if our 
%   bearing is split into two bearings by a ring groove)
% - n_x_noTaylor = circumferential number of nodes for the FVM model 
%   (method=1) or the SBFEM model with a call of the eigensolver at every
%   time step (method=2); the SBFEM model with Taylor approximations of the
%   eigenvalues and eigenvectors (method=3) ignores this input variable
%   because its number of nodes was already defined in 'A_precomputation.m'
%   [-]

% OUTPUT VARIABLES
% - t = column vector containing the time values at to the output steps [s]
% - z = the i-th row of this matrix describes the state-space vector at
%   the i-th output step; the state space vector includes the horizontal 
%   displacement [m], vertical displacement [m], horizontal velocity [m/s], 
%   and vertical velocity [m/s] of the shaft (in this order)
% - n_calls = how many times the Reynolds equation was solved [-]



% define grid


if method == 3                                                             % if the SBFEM with Taylor approximations of the eigenvalues and eigevectors is used
    load('3_coefficients.mat','n_x')                                       % the circumferential number of nodes n_x was already defined in 'A_precomputation.m'
else                                                                       % if a different simulation method is used
    n_x = n_x_noTaylor;                                                    % the circumferential number of nodes n_x was defined in 'B_computation.m' and stored in the variable n_x_noTaylor
end
n_y = round((l_b/2)/(pi*d_b)*n_x)+1;                                       % number of nodes in the axial direction (for one half of the bearing, i.e., l_b/2) in the FVM solution; these are also used as integration points for the bearing forces in the SBFEM, unless this integration is performed analytically



% load eigenvalue and eigenvector derivatives if required


if method == 3                                                             % if the SBFEM model with eigenvalue and eigenvector derivatives is used
    
    load('1_locations.mat','constr_vec','switch_vec')                      % constr_vec contains the points throughout the range of relative eccentricities where Taylor series' will be constructed and ...
    load('2_reductions.mat','red_vec')                                     % ... red_vec describes, for each of these points, the size of the subset of modes to be considered in the computation of the pressure field ("red" for reduction, as in modal reduction)
    load('3_coefficients.mat','n_tay','gamma_ref',...                      % load all data necessary for the Taylor series approximations of the eigenvalues and eigenvectors (this data was generated by 'A_precomputation.m')
        'Vd_allpoints_mat','Ld_allpoints_vec','ind_vec')
    Ld_allpoints_vec = Ld_allpoints_vec * ...                              % adjust eigenvalues and their derivatives to the defined bearing dimensions
        ((l_b/d_b)/gamma_ref)^2;
    
else
    
    Vd_mat = [];                                                           % no eigenvector derivatives are needed, define as empty variable
    Ld_vec = [];                                                           % no eigenvalue derivatives are needed, define as empty variable
    epsilon_constr = [];                                                   % no Taylor series is constructed, define point of construction as empty variable
    n_tay = [];                                                            % no Taylor series is constructed, define oder of Taylor series as empty variable
    
end



% derivation of some parameters


Delta_t = t_end-t_start;                                                   % time span of run-up
Delta_omega = omegaspan(1,2)-omegaspan(1,1);                               % difference in angular velocity between beginning and end of run-up
omega_dot = Delta_omega/Delta_t;                                           % angular acceleration of shaft

M_inv = diag([1,1,m,m].^(-1));                                             % inverse mass matrix
K = [0,0,-1,0;0,0,0,-1;0,0,0,0;0,0,0,0];                                   % matrix which couples the velocities in the state-space vector to their counterparts in the derivative of the state-space vector

n_out = fac_frq*(t_end-t_start)*max(abs(omegaspan))/(2*pi);                % number of output steps
tspan = linspace(t_start,t_end,n_out);                                     % vector containing the time values corresponding to the output steps

n_calls = 0;                                                               % this variable will count the number of calls of the Reynolds equation, it is initialized as zero



% time integration


options = odeset('RelTol',1e-9);
[t,z] = ode23t(@(t,z) odefun(t,z),tspan,z_0,options);                      % call Matlab's ode23t time integrator to solve equation of motion in the time domain



% define function describing the equation of motion


function dzdt = odefun(t,z)                                                % this function describes the equation of motion; it computes the forces acting on the shaft and then computes the derivative of the state-space vector
    
    % analyze kinematic variables of the shaft
    
    x = z(1,1);                                                            % horizontal shaft position is extracted from state-space vector
    y = z(2,1);                                                            % vertical shaft position is extracted from state-space vector
    x_dot = z(3,1);                                                        % horizontal shaft velocity is extracted from state-space vector
    y_dot = z(4,1);                                                        % vertical shaft velocity is extracted from state-space vector
    if x == 0 && y == 0                                                    % if the eccentricity is zero
        x = eps;                                                           % set shaft displacements to tiny nonzero value to avoid singularity in the computations below
        y = eps;                                                           % set shaft displacements to tiny nonzero value to avoid singularity in the computations below
    end
    q = sqrt(x^2+y^2);                                                     % absolute eccentricity
    q_dot = (y_dot*y+x_dot*x)/q;                                           % rate of change of absolute eccentricity
    X_att = atan2(y,x);                                                    % attitude angle
    X_att_dot = (y_dot*x-x_dot*y)/q^2;                                     % rate of change of attitude angle
    epsilon = q/c;                                                         % relative eccentricity
    epsilon_dot = q_dot/c;                                                 % rate of change of relative eccentricity
    omega = omegaspan(1,1)+omega_dot*(t-t_start);                          % angular velocity of shaft
    phi = phi_u+omegaspan(1,1)*(t-t_start)+0.5*omega_dot*(t-t_start)^2;    % angular position of unbalance
    
    % solve Reynolds equation, compute hydrodynamic forces
    
    if method == 1                                                         % if the FVM is chosen
        [f_b_x,f_b_y] = fvm(d_b,l_b,c,omega,mu,n_x,n_y,epsilon,...         % solve Reynolds equation with FVM, store components of bearing force in f_b_x and f_b_y
            epsilon_dot,X_att,X_att_dot,pb1,pb2);
    elseif method == 2                                                     % if the SBFEM with solution of eigenvalue problem chosen
        red = n_x;                                                         % no modal reduction
        flag_Tay = 0;                                                      % disable Taylor series
        [f_b_x,f_b_y] = sbfem(d_b,l_b,c,omega,mu,n_x,n_y,epsilon,...       % solve Reynolds equation with SBFEM, store components of bearing force in f_b_x and f_b_y   
            epsilon_dot,X_att,X_att_dot,epsilon_constr,n_tay,...
            Vd_mat,Ld_vec,flag_Tay,red,pb1,pb2);
    elseif method == 3                                                     % if the SBFEM with eigenvalue and eienvector derivatives is chosen 
        [~,k] = min(abs(switch_vec-epsilon));
        if epsilon < switch_vec(k,1)
            k = k-1;
        end
        epsilon_constr = constr_vec(k,1);                                  % point where Taylor series was constructed
        red = red_vec(k,1);                                                % number of modes to consider  
        ind = ind_vec(k,1);
        Vd_mat = Vd_allpoints_mat(:,ind:(ind-1+red*(n_tay+1)));            % eigenvector derivatives
        Ld_vec = Ld_allpoints_vec(1,ind:(ind-1+red*(n_tay+1)));            % eigenvalue derivatives
        flag_Tay = 1;                                                      % enable Taylor series
        [f_b_x,f_b_y] = sbfem(d_b,l_b,c,omega,mu,n_x,n_y,epsilon,...       % solve Reynolds equation with SBFEM, store components of bearing force in f_b_x and f_b_y
            epsilon_dot,X_att,X_att_dot,epsilon_constr,n_tay,...
            Vd_mat,Ld_vec,flag_Tay,red,pb1,pb2);
    end
    
    % summation of forces, computation of derivative of state-space vector
    
    f_shaft_x = -f_b_x*mbf+omega^2*unb*cos(phi);                           % horizontal force acting on the shaft, considering bearing force and unbalance
    f_shaft_y = -f_b_y*mbf+omega^2*unb*sin(phi)-m*g;                       % vertical force acting on the shaft, considering bearing force, unbalance, and gravitation
    rhs = [0;0;f_shaft_x;f_shaft_y];                                       % right-hand side vector of equation of motion
    dzdt = M_inv*(-K*z+rhs);                                               % equation of motion is solved for the derivative of the state-space vector
    n_calls = n_calls+1;                                                   % count number of calls of the Reynolds equation
    
end



end