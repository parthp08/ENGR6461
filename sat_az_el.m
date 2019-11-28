function [az, el] = sat_az_el(P_sat,P_u)
    
    xs = P_sat(1);
    ys = P_sat(2);
    zs = P_sat(3);
    
    xu = P_u(1);
    yu = P_u(2);
    zu = P_u(3);
    
	[lat, lon, h] = ECEF2WGS(P_u, 0);
    
    % vector from user to sat
    X = [xs-xu; ys-yu; zs-zu];

    % Transformation from the ECEF to ENU
    lon_s = sin(lon);
    lon_c = cos(lon);
    lat_s = sin(lat);
    lat_c = cos(lat);
    R_L = [
        -lon_s, lon_c, 0;
        -lat_s*lon_c, -lat_s*lon_s, lat_c;
        lat_c*lon_c, lat_c*lon_s, lat_s
    ];

    % X_L % representation of X in ENU
    X_L = R_L*X;
    x_e = X_L(1); % east position
    x_n = X_L(2); % north position
    x_u = X_L(3); % up position

    az = atan2(x_e, x_n);

    el = asin(x_u / sqrt(x_e^2 + x_n^2 + x_u^2));

end
