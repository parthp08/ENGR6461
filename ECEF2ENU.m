%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Converts vector ocordinates from ECEF frame to WGS84 frame
%
%%% references
% -------------
%[1] 'Global Positioning System: Signals, Measurements, and 
%     Performance', Pratap Misra, Per Enge.
%
%%% inputs
% ----------
% P : array, size(3,1), vector in ECEF Frame, meters
%
%%% outputs
%-----------
% X_L : array, size(3,1), vector in ENU Frame, meters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X_L = ECEF2ENU(P)
    
    [lat, lon, ~] = ECEF2WGS(P, 0);
    
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
    X_L = R_L*P;
    end