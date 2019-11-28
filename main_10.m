clc;
load project_data.mat;

% with all errors

% using the answer of main_p as an initial position
Pu_0 = [3.5060e+06; 0.3280e6; 4.8906e6]; % initial user position % origin
cb_u0 = -3.2134e+05; % initial estimation of user clock bias
t_rcv = iono(1); % receiving GPS time of the week
tau = 0.075; % intial estimation of transmission time between sat and user
c = 299792458;   % m/sec % speed of light

% satellites' position(P) (6 satellites => a to f)
P_a = sat_position(eph(:,1), t_rcv, tau); % can take these out of the loop if dont want to recalculate tau
P_b = sat_position(eph(:,2), t_rcv, tau);
P_c = sat_position(eph(:,3), t_rcv, tau);
P_d = sat_position(eph(:,4), t_rcv, tau);
P_e = sat_position(eph(:,5), t_rcv, tau);
P_f = sat_position(eph(:,6), t_rcv, tau);

P_sat_arr = [P_a, P_b, P_c, P_d, P_e, P_f];


%%% corrections for pseudorange (for all satellites)
% satellite clock corrections
delta_t_sv = [];
delta_t_sv(1) = sat_clock_error(eph(:,1), t_rcv, tau);
delta_t_sv(2) = sat_clock_error(eph(:,2), t_rcv, tau);
delta_t_sv(3) = sat_clock_error(eph(:,3), t_rcv, tau);
delta_t_sv(4) = sat_clock_error(eph(:,4), t_rcv, tau);
delta_t_sv(5) = sat_clock_error(eph(:,5), t_rcv, tau);
delta_t_sv(6) = sat_clock_error(eph(:,6), t_rcv, tau);

% ionosphere and troposphere errors
T_amb = 15;  % deg C
P_amb = 101.325; % kPa
P_vap = 0.85;    % kPa
[az1, el1] = sat_az_el(P_sat_arr(:,1), Pu_0);
[az2, el2] = sat_az_el(P_sat_arr(:,2), Pu_0);
[az3, el3] = sat_az_el(P_sat_arr(:,3), Pu_0);
[az4, el4] = sat_az_el(P_sat_arr(:,4), Pu_0);
[az5, el5] = sat_az_el(P_sat_arr(:,5), Pu_0);
[az6, el6] = sat_az_el(P_sat_arr(:,6), Pu_0);
[lat_u, lon_u, h] = ECEF2WGS(Pu_0, 0);
I_e = [];
I_e(1) = iono_error(iono, lat_u, lon_u, el1, az1);
I_e(2) = iono_error(iono, lat_u, lon_u, el2, az2);
I_e(3) = iono_error(iono, lat_u, lon_u, el3, az3);
I_e(4) = iono_error(iono, lat_u, lon_u, el4, az4);
I_e(5) = iono_error(iono, lat_u, lon_u, el5, az5);
I_e(6) = iono_error(iono, lat_u, lon_u, el6, az6);
T_e = [];
T_e(1) = tropo_error(T_amb, P_amb, P_vap, el1);
T_e(2) = tropo_error(T_amb, P_amb, P_vap, el2);
T_e(3) = tropo_error(T_amb, P_amb, P_vap, el3);
T_e(4) = tropo_error(T_amb, P_amb, P_vap, el4);
T_e(5) = tropo_error(T_amb, P_amb, P_vap, el5);
T_e(6) = tropo_error(T_amb, P_amb, P_vap, el6);

% apply correction    
pr_corrected = pr' + (c.*(delta_t_sv'))- cb_u0 - I_e' - T_e';
pr = pr_corrected';

%%% Least Square Solution
[P_u, cb_u, delta_r] = least_square_sol(P_sat_arr, Pu_0, pr);

% assign new value for initial position
Pu_0 = P_u;
%cb_u0 = cb_u;

% On 3D Earth
pos_x = [P_u(1), P_a(1), P_b(1), P_c(1), P_d(1), P_e(1), P_f(1)];
pos_y = [P_u(1), P_a(2), P_b(2), P_c(2), P_d(2), P_e(2), P_f(2)];
pos_z = [P_u(1), P_a(3), P_b(3), P_c(3), P_d(3), P_e(3), P_f(3)];

earth_sphere('m')
hold on
plot3(pos_x, pos_y, pos_z, 'o')
