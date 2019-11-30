%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute the ionospheric errors in pseudorange mesurement
%
%%% references
% -------------
%[1] https://moodle.concordia.ca/moodle/pluginfile.php/3799910/mod_resource/content/1/Project_Appendix.pdf
%[2] https://moodle.concordia.ca/moodle/pluginfile.php/3799909/mod_resource/content/1/ENGR6461_project.pdf
%[3] 'Global Positioning System: Theory and Applications'; edited by 
%     Parkinson, Spilker, Axelrad, Enge; AIAA; 1996
%[4] 'Aerospace Navigation Systems'; edited by A.V.Nebylv, J.Watson
%
%%% inputs
% ----------
% iono_data : array, size(9,1), data with GPS time for measurement and 
%             eight parameters for ionosphere delay estimation
% lat_u : float, latitude of user, rad
% lon_u : float, longitude of user, rad
% E_s : float, elevation of sat from the user local tangent plane, rad
% A_s : float, azimuth of sat from the user local tangent plane, rad
%
%%% outputs
%-----------
% I_e : float, ionospheric delay, meters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I_e = iono_error(iono_data, lat_u, lon_u, E_s, A_s)

    pi = 3.1415926535898;    % GPS standard value for pi

    phi_u = lat_u;
    lamda_u = lon_u;

    % Ref [2] page 2
    tow = iono_data(1);
    a0 = iono_data(2);
    a1 = iono_data(3);
    a2 = iono_data(4);
    a3 = iono_data(5);
    b0 = iono_data(6);
    b1 = iono_data(7);
    b2 = iono_data(8);
    b3 = iono_data(9);

    %%% Klobuchar model for ionosphere  % Ref [1] page 4-6 and Ref [2] page 146-147 and Ref [4] page 78-79
    % Earth-centered angle between user and IP(ionosphere point % Ref[4] page 78)
    psi_ip = (0.0137/(0.11 + E_s)) - 0.022; % semicircles
    psi_ip = pi * psi_ip; % to convert to rad
    
    % subionospheric latitide
    phi_ip = phi_u + psi_ip*cos(A_s);    % semicircles
    phi_ip = pi * phi_ip; % to convert to rad
    if phi_ip > 0.416
        phi_ip = 0.416;
    end
    if phi_ip < -0.416
        phi_ip = -0.416;
    end
    
    % subionospheric longitude
    lamda_ip = lamda_u + psi_ip*sin(A_s)/cos(phi_ip);  % semicircles
    lamda_ip = pi * lamda_ip; % to convert to rad

    % geomagnetic lattitude of the subionospheric location looking toward each GPS satellite
    phi_m = phi_ip + 0.064*cos(lamda_ip-1.617);   % semicircles
    phi_m = pi * phi_m; % to convert to rad

    % local time(t) at the subionospheric point
    t_gps = mod(tow, 86400);  % Appendix % moodle % second last page  % written note
    t_ip = 4.32e4*lamda_ip + t_gps; % seconds
    if t_ip > 86400
        t_ip = t_ip - 86400;
    end
    if t_ip < 0
        t_ip = t_ip + 86400;
    end
    
    % obliquity factor / slant factor (dimensionless)
    % to convert to slant time delay 
    F_ip = 1.0 + 16.0*((0.53-E_s)^3);

    % Amplitude % in seconds
    AM_ip = a0 + a1*(phi_m^1) + a2*(phi_m^2) + a3*(phi_m^3);
    if AM_ip < 0
        AM_ip = 0;
    end

    % Period % in seconds
    PER_ip = b0 + b1*(phi_m^1) + b2*(phi_m^2) + b3*(phi_m^3);
    if PER_ip < 72000
        PER_ip = 72000;
    end

    % Phase % in rad
    X_ip = 2*pi*(t_ip - 50400) / PER_ip;

    % Ionospheric Time delay (T_iono)  % in seconds
    if abs(X_ip) >= 1.57
        T_iono = F_ip*5e-9;
    else   % for abs(X_ip) < 1.57
        T_iono = F_ip*(5e-9 + AM_ip*(1 - (X_ip^2)/2 + (X_ip^4)/24));
    end

    % Ionosphere Error % in meters
    c = 299792458;   % speed of light(GPS constant) % m/s
    I_e = c * T_iono;    % m

end

