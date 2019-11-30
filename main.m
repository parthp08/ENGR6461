%%% main script for calculating the user position from GPS data
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% user position calculation
%%% without considering any errors in the pseudoranges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPS data recieved from the satellite
load project_data.mat;

% obtain receiving GPS time of the week from the data 
t_rcv = iono(1);
tau = 0.075; % intial estimation of transmission time between sat and user
c = 299792458; % m/sec % constant % speed of light

% Calculate the satellites position from the ephemeris data
P_a = sat_position(eph(:,1), t_rcv, tau);
P_b = sat_position(eph(:,2), t_rcv, tau);
P_c = sat_position(eph(:,3), t_rcv, tau);
P_d = sat_position(eph(:,4), t_rcv, tau);
P_e = sat_position(eph(:,5), t_rcv, tau);
P_f = sat_position(eph(:,6), t_rcv, tau);

P_sat_arr = [P_a, P_b, P_c, P_d, P_e, P_f]; % sat pos array

% inital user position and clock bias assumption
% at the centre of the earth in ECEF frame
Pu_0 = [0; 0; 0];
bu_0 = 0;

% Calculate the Inital Good Estimate of the inital Position
% to use later for error correction

% just to get loop started
delta_r = [100;100;100];
while norm(delta_r(1:3)) > 1e-5

    % Least Square Solution
    delta_r = least_square_sol(P_sat_arr, Pu_0, pr);

    % calulate user position and user clock bias
    Pu_0 = delta_r(1:3) + Pu_0;
    bu_0 = delta_r(4);
end

% user position and clock error % without any errors considered
Pu = delta_r(1:3) + Pu_0;
bu = delta_r(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% user position calculation
%%% with all available errors considered
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use improved initial estimation for incorporating errors into the 
% mesurements
Pu_0 = Pu;
bu_0 = bu;

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
[lat_u, lon_u, ~] = ECEF2WGS(Pu_0, 0);
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
pr_corrected = pr' + (c.*(delta_t_sv'))- bu_0 - I_e' - T_e';
pr = pr_corrected';

%%% Least Square Solution
delta_r = least_square_sol(P_sat_arr, Pu_0, pr);

% user position and clock error % with errors considered
Pu = delta_r(1:3) + Pu_0;
bu = delta_r(4);

% user position in Longitude and Lattitude
[lat, lon, h] = ECEF2WGS(Pu, 1);  % 1 = units in deg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the user and satellites positions on ECEF frame
sat_pos_x = [P_a(1), P_b(1), P_c(1), P_d(1), P_e(1), P_f(1)];
sat_pos_y = [P_a(2), P_b(2), P_c(2), P_d(2), P_e(2), P_f(2)];
sat_pos_z = [P_a(3), P_b(3), P_c(3), P_d(3), P_e(3), P_f(3)];
figure();
earth_sphere('m'); % earth 3d map
hold on;
scatter3(sat_pos_x, sat_pos_y, sat_pos_z, 'o')
scatter3(Pu(1), Pu(2), Pu(3), '*', 'w');

figure();
geoscatter(lat,lon, '*', 'r'); % on 2D map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
