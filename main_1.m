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

% just to get loop started
delta_r = [100; 100; 100; 100];
while norm(delta_r(1:3)) > 1e-5 % delta_r convergence to zero (user location == initial location)

    %%% corrections for pseudorange (for all satellites)
    % satellite clock corrections
    delta_t_sv = [];
    for i = 1:6 % 6 = number of satellites
        delta_t_sv(i) = sat_clock_error(eph(:,i), t_rcv, tau);
    end
    % ionosphere and troposphere errors
    T_amb = 15;  % deg C
    P_amb = 101.325; % kPa
    P_vap = 0.85;    % kPa
    I_e = [];
    T_e = [];
    for i = 1:6
        [az, el] = sat_az_el(P_sat_arr(:,i), Pu_0);
        [lat_u, lon_u, h] = ECEF2WGS(Pu_0, 0);
        I_e(i) = iono_error(iono, lat_u, lon_u, el, az);
        T_e(i) = tropo_error(T_amb, P_amb, P_vap, el);
    end
    % apply correction    
    pr_corrected = pr' + (c.*(delta_t_sv')) - I_e' - T_e';
    pr = pr_corrected';
    
    %%% Least Square Solution
    [P_u, cb_u, delta_r] = least_square_sol(P_sat_arr, Pu_0, pr);
    
    % assign new value for initial position
    Pu_0 = P_u;
    %cb_u0 = cb_u;
    disp("parth");  % just to see how many times loop runs :) :)
end

% somethings wrong 
