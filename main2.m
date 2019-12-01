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
P_sat_arr = sat_position(eph, t_rcv, tau); % sat pos array

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
delta_t_sv = sat_clock_error(eph, t_rcv, tau);

% ionosphere and troposphere errors
T_amb = 15;  % deg C
P_amb = 101.325; % kPa
P_vap = 0.85;    % kPa
[az, el] = sat_az_el(P_sat_arr, Pu_0);
[lat_u, lon_u, ~] = ECEF2WGS(Pu_0, 0);
I_e = iono_error(iono, lat_u, lon_u, el, az);
T_e = tropo_error(T_amb, P_amb, P_vap, el);

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
figure();
earth_sphere('m'); % earth 3d map
hold on;
scatter3(P_sat_arr(1,:),P_sat_arr(2,:),P_sat_arr(3,:), 'o')
scatter3(Pu(1), Pu(2), Pu(3), '*', 'w');

figure();
geoscatter(lat,lon, '*', 'r'); % on 2D map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
