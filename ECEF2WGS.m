function [phi, lambda, h] = ECEF2WGS(P, deg)
 
    x = P(1);
    y = P(2);
    z = P(3);
    % WGS84 Model
    a = 6378137.0;  % meters % semi-major axis
    b = 6356752.314;  % meters % semi-minor axis
    e = sqrt(1 - (b^2 / a^2));% eccentricity
    
    %%% Step 1
    if (x^2 + y^2) == 0  % if True algorithm stops here
        if z > 0
            phi = pi/2;
            h = z - b;
        else % z < 0
            phi = -pi/2;
            h = -z - b;
        end
        % return "Longitude cannot be defined", phi, h
        lambda = 0;
     
    else
        
    %%% Step 2
    % assuming x^2 + y^2 > 0
    if x == 0
        lambda_ = pi/2;
    elseif x > 0
        lambda_ = atan(y/x);
    elseif x < 0 && y >= 0 
        lambda_ = pi + atan(y/x);
    elseif x < 0 && y < 0
        lambda_ = -pi + atan(y/x);
    end
        
    %%% Step 3
    % sensitivity levels
    e_phi = 0.01;
    e_h = 0.01;

    % starting values
    phi_s = pi/2;
    d_s = a / sqrt(1 - e^2);
    h_s = 0;

    %%% Step 4
    % computing new values related to the latitude
    phi_n = atan2(z*(d_s+h_s) , ((sqrt(x^2 + y^2))*((d_s*(1 - e^2))+h_s)));
    d_n = a / sqrt(1 - ((e*sin(phi_n))^2));

    %%% Step 5
    % computing new altitude value
    if phi_n <= pi/4 || phi_n >= -pi/4
        h_n = (sqrt(x^2 + y^2)/cos(phi_n)) - d_n;
    end
    if phi_n > pi/4 || phi_n < -pi/4
        h_n = (z/sin(phi_n)) - (d_n*(1 - e^2));
    end
    
    %%% Step 6
    % Evaluating the convergence of the obtained value for phi
    while abs(phi_n - phi_s) >= e_phi && abs(h_n - h_s) >= e_h % if not converged
        phi_s = phi_n;
        d_s = d_n;
        h_s = h_n;

        %%% Step 4
        % computing new values related to the latitude
        phi_n = atan2(z*(d_s+h_s), (sqrt(x^2 + y^2))*(d_s*(1 - e^2) + h_s));
        d_n = a / sqrt(1 - ((e*sin(phi_n))^2));

        %%% Step 5
        % computing new altitude value
        if phi_n <= pi/4 || phi_n >= -pi/4
            h_n = (sqrt(x^2 + y^2)/cos(phi_n)) - d_n;
        end
        if phi_n > pi/4 || phi_n < -pi/4
            h_n = (z/sin(phi_n)) - (d_n*(1 - e^2));
        end
    end

    % when converged
    % abs(phi_n - phi_s) >= e_phi and abs(h_n - h_s) >= e_h
    lambda = lambda_;   % in radians
    phi = phi_n;
    h = h_n;
    
    if deg == 1
        lambda = rad2deg(lambda_);   % in degrees
        phi = rad2deg(phi_n);
        h = rad2deg(h_n);
    end
    
    end
end

