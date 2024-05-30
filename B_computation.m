% AUTHOR, AFFILIATION, DATE
% Simon Pfeil, OvGU Magdeburg, 30.05.24

% DESCRIPTION
% Run-up simulation with a rotor represented by a point mass in a plain 
% bearing; the Reynolds equation is solved with either the FVM or the 
% SBFEM; this is the main script where the parameters are defined and 
% which then calls the function 'run_up.m' to perform the simulation;
% run this script in Matlab



% clear variables, clear console, close figures, etc. ...


clear variables
close all
clc
addpath(genpath(pwd()))
dbstop if error



% computational parameters that were not defined in 'A_precomputation.m'


method = 3;                                                                % computational method for Reynolds equation: (1) = FVM, (2) = SBFEM with solution of eigenvalue problem at every time step, (3) = SBFEM with eigenvalue derivatives []
n_x_noTaylor = 100;                                                        % circumferential number of nodes (will be ignored if method=3, because in that case, n_x was already defined in 'A_precomputation.m') []
fac_frq = 20;                                                              % frequency for saving the results relative to the max frequency of rotation []
mbf = 1;                                                                   % multiplier for the bearing forces; use 1 for a plain bearing or 2 for a bearing with circumferential groove (which is interpreted as two separate but identical bearings) []



% physical/technical parameters


d_b = 0.1;                                                                 % bearing diameter [m]
l_b = 0.1;                                                                 % bearing length (for bearings with ring groove, the two segments separated by the groove are treated as individual bearings; in that case, set l_b to the length of such a segment) [m]
mu = 0.02;                                                                 % viscosity [Pa*s]
unb = 0.005;                                                               % unbalance [kg*m]
phi_u = 0;                                                                 % angular position of the unbalance, relative to the rotation of the shaft [rad]
g = 9.81;                                                                  % gravitation [N/kg]
t_start = 0;                                                               % time at the beginning [s]
z_0 = [0;0;0;0];                                                           % initial conditions of the rotor (horizontal displacement, vertical displacement, horizontal velocity, vertical velocity), i.e., state-space vector at t=t_start [m]/[m/s]
pb1 = 0;                                                                   % pressure at the first bearing boundary (zero is atmosphere) [Pa]
pb2 = 0;                                                                   % pressure at the second bearing boundary (zero is atmosphere) [Pa]
c = 150e-6;                                                                % radial clearance [m]
m = 200;                                                                   % rotor mass (rotor=shaft) [kg]
omegaspan = 2*pi*[60,140];                                                 % angular velocities at the beginning and at the end, assuming constant acceleration [rad/s]
t_end = 1;                                                                 % time at the end [s]



% execute simulation


tic0 = tic;                                                                % start clock for computational time

[t,z,n_calls] = run_up(d_b,l_b,c,mu,m,omegaspan,t_start,t_end,...          % execute simulation; the time values and state-space vectors are stored in t and z, respectively, for all output steps; the number of calls of the Reynolds equation is stored in n_calls
    z_0,method,g,unb,fac_frq,phi_u,pb1,pb2,mbf,n_x_noTaylor);

t_comp = toc(tic0);                                                        % stop clock for compuational time



% store and display results


disp(['Elapsed time is ',num2str(t_comp),' seconds.'])
disp(['Reynolds equation was solved ',num2str(n_calls),' times.'])
disp(['Maximum relative eccentricity is ',num2str(max(sqrt(z(:,1).^2+...
    z(:,2).^2))/c)])

n_out = length(t);                                                         % number of output steps

% save(['workspace_method_',num2str(method)])                                % save workspace

figure(1)
plot(t(:,1),z(:,2)/c)                                                      % plot vertical shaft oscillation over time
hold on
