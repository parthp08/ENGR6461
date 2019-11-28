clc;
load project_data.mat;

Pu_0 = [0; 0; 0]; % initial user position % origin
t_rcv = iono(1); % receiving GPS time of the week
tau = 0.075; % intial estimation of transmission time between sat and user

% satellites' position(P) (6 satellites => a to f)
P_a = sat_position(eph(:,1), t_rcv, tau); % can take these out of the loop if dont want to recalculate tau
P_b = sat_position(eph(:,2), t_rcv, tau);
P_c = sat_position(eph(:,3), t_rcv, tau);
P_d = sat_position(eph(:,4), t_rcv, tau);
P_e = sat_position(eph(:,5), t_rcv, tau);
P_f = sat_position(eph(:,6), t_rcv, tau);

P_sat_arr = [P_a, P_b, P_c, P_d, P_e, P_f];

% just to get loop started
delta_r = [100; 100; 100; 100];
while norm(delta_r(1:3)) > 1e-5 % delta_r convergence to zero (user location == initial location)

    % Least Square Solution
    [P_u, cb_u, delta_r] = least_square_sol(P_sat_arr, Pu_0, pr);

    % assign new value for initial position
    Pu_0 = P_u;
    %disp("parth");  % just to see how many times loop runs :) :)
end

% Most Promising result till now (Just see in plots with Earth Model)
% Lat = 54.438  % deg
% Lon = 5.34469 % deg
% h = -337619.9 % m

