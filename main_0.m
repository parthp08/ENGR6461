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

% unit vectors(e) from user position to satellites' position
% https://www.mathworks.com/matlabcentral/answers/45376-create-cartesian-unit-vectors-from-2-points
e_a = unit_vector(Pu_0, P_a);
e_b = unit_vector(Pu_0, P_b);
e_c = unit_vector(Pu_0, P_c);
e_d = unit_vector(Pu_0, P_d);
e_e = unit_vector(Pu_0, P_e);
e_f = unit_vector(Pu_0, P_f);

% form H matrix
H = [-e_a', 1;...
     -e_b', 1;...
     -e_c', 1;...
     -e_d', 1;...
     -e_e', 1;...
     -e_f', 1];

H_terms = (inv(H'*H))*H';

% pseudorange %delta_rho = rho - rho_0
% rho => pseudorange of user from satellite (given in data)
% rho_0 => pseudorange of user initial location from satellite
% https://www.mathworks.com/matlabcentral/answers/16848-how-to-find-distance-between-two-points
rho = pr';
rho_0 = [norm(P_a - Pu_0);
         norm(P_b - Pu_0);
         norm(P_c - Pu_0);
         norm(P_d - Pu_0);
         norm(P_e - Pu_0);
         norm(P_f - Pu_0)];
delta_rho = rho - rho_0;

% Least Square Solution for user position(P_u)
% delta_r = P_u - Pu_0
delta_r = H_terms*delta_rho;
P_u = delta_r(1:3) + Pu_0;

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

% with initial user position (52.3667° N, 4.8945° E) and h = 6371000
% that's in ECEF = (7764.478e3, 664.9e3, 10073.238e3)
% xu = 4.271071879434444e+06
% yu = 4.365719402779287e+05
% zu = 5.691189381888778e+06
% Lat, Lon = 53.13, 5.84 (in Netherland => De Alde Feanen National Park)
% h = 764515.2

% with initial user position (from above answer)
% (xu0,yu0,zu0) = (4.2711e6; 0.4366e6; 5.6912e6)
% xu = 3.532776465952689e+06
% yu = 3.324967171598776e+05
% zu = 4.916607519499978e+06
% Lat, Lon = 54.37, 5.38 (back to the sea :(  )
% h = -300702.5



