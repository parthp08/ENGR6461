% Calculating the position of the user with Least square solution
% without considering any errors

% assuming initial conditions to zero xu0,yu0,zu0 = (0;0;0) % origin
% rho_0 = from origin to sat_pos
% rho = from data
% delta_rho = rho - rho_0
% delta_r = r - r_o  == r   % r_0 == origin
% r => [xu;y;,zu;bcu]
% H-matrix elements => unit vectors % sat_pos from origin

clc;
load project_data.mat;

% so that we are getting location in sea near Netherland, what if we use 
% location near it to see if we are getting any better results than before.
%Pu_0 = [7764.478e3; 664.9e3; 10073.238e3]; % Amsterdam location
%Pu_0 = [4.2711e6; 0.4366e6; 5.6912e6]; % answer of above initial location
Pu_0 = [0; 0; 0]; % initial user position % origin
t_rcv = iono(1); % receiving GPS time of the week
tau = 0.075; % intial estimation of transmission time between sat and user

% satellites' position(P) (6 satellites => a to f)
P_a = sat_position(eph(:,1), t_rcv, tau);
P_b = sat_position(eph(:,2), t_rcv, tau);
P_c = sat_position(eph(:,3), t_rcv, tau);
P_d = sat_position(eph(:,4), t_rcv, tau);
P_e = sat_position(eph(:,5), t_rcv, tau);
P_f = sat_position(eph(:,6), t_rcv, tau);

P_sat_arr = [P_a, P_b, P_c, P_d, P_e, P_f];

% Least Square Solution
[P_u, cb_u, delta_r] = least_square_sol(P_sat_arr, Pu_0, pr);

% to convert ECEF to LLH or Visa Versa
% https://www.oc.nps.edu/oc2902w/coord/llhxyz.htm
% x,y,z are in Km 
% Lat, Lon in deg and h in meters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% result analysis
%%%%%%%%%%%%%%%%%%%%%%%%%

% with initial user position (0, 0, 0)
% xu = 4.287476159386505e+06
% yu = 4.276037784550767e+05
% zu = 5.905701167290296e+06
% lat, lon = 54.042, 5.69542 (in sea near Netherland)
% h = 946279.1



