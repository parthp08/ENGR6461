function P = sat_position(eph_data, t_rcv, tau)

    %Calculate the Position of Satellite in ECEF frame given the ephemeris data of the satellite
    %and GPS receiving time of the week

    %Inputs
    %-------
    %eph_data : array, shape(21,1), ephemeris data for a satellite
    %t_rcv : float, receiving GPS time of the week, seconds
    %tau: float, transmission time between sat and user, seconds

    %Returns
    %--------
    %P = [x; y; z] : array, floats, Position of Satellte in ECEF frame, in meters
    
    %References
    %-----------
    %[1] https://moodle.concordia.ca/moodle/pluginfile.php/3799910/mod_resource/content/1/Project_Appendix.pdf
    %[2] https://moodle.concordia.ca/moodle/pluginfile.php/3799909/mod_resource/content/1/ENGR6461_project.pdf
    %[3] 'Global Positioning System: Theory and Applications'; edited by Parkinson, 
    %     Spilker, Axelrad, Enge; AIAA; 1996
    %[4] 'Aerospace Navigation Systems'; edited by A.V.Nebylv, J.Watson
    %[5] https://ascelibrary.org/doi/pdf/10.1061/9780784411506.ap03
    %[6] GNSS Applications and Methods; Scott Gleason, Demoz Gebre-Egzibher
    
    pi = 3.1415926535898;    % GPS standard value for pi
    
    % eph parameters (Ref [2], page 2)
    cr_sin = eph_data(1);
    delta_n = eph_data(2) * pi;
    M_0 = eph_data(3) * pi;
    cu_cos = eph_data(4);
    e = eph_data(5);
    cu_sin = eph_data(6);
    sqrta = eph_data(7);
    t_0e = eph_data(8);
    ci_cos = eph_data(9);
    omg0 = eph_data(10) * pi;
    ci_sin = eph_data(11);
    i_0 = eph_data(12) * pi;
    cr_cos = eph_data(13);
    w = eph_data(14) * pi;
    omg_dot = eph_data(15) * pi;
    i_dot = eph_data(16) * pi;
    %tgd = eph_data(17);
    %toc = eph_data(18);
    %af2 = eph_data(19);
    %af1 = eph_data(20);
    %af0 = eph_data(21);
    
    %%% Satellite Position Calculation Algorithm (Ref[1] page 3, Ref [3] page 138)
    % constants
    mu = 3.986005e14;    % m^3 / s^2 % Earth universal gravitational parameter
    omg_e_dot = 7.2921151467e-5; % Earth Rotation Rate
    A = sqrta^2;    % semimajor axis
    
    n0 = sqrt(mu / A^3); % computed mean motion % rad/sec
    
    %tau = 0.075;
    t = t_rcv - tau; % t == t_tr
    t_k = t - (t_0e);  % time from ephemeris reference time(t_0e)
    % from the footnote(page 138 Ref [3])
    % accounting for beginning or end of week crossovers
    if t_k > 302400  % seconds
        t_k = t_k - 604800;
    end
    if t_k < -302400
        t_k = t_k + 604800;
    end
        
    n = n0 + delta_n;    % corrected mean motion
    
    M_k = M_0 + n*t_k;   % Mean anomaly
        
    % (Ref [1] page 3(Note))
    % E_k = M_k + e*np.sin(E_k)    % Kepler's Equ for the eccentric anomaly E_k % rad
    E_old = M_k;
    E_new = M_k + e*sin(E_old);
    while E_new - E_old >= 1e-8  % stop iteration when E_new - E_old is less than 1e-8
        E_old = E_new;
        E_new = M_k + e*sin(E_old);
    end
    E_k = E_new; % Essentric anomaly at time t_k
    
    % True anomaly as a function of E_k at time t_k
    v_k = atan2(sqrt(1-e^2)*sin(E_k), cos(E_k)-e);

    E_k = acos((e+cos(v_k)) / (1 + e*cos(v_k)));    % Eccentric anomaly

    phi_k = v_k + w;   % Argument of latitude
    
    % Second harmonic perturbations
    phi_k_2 = 2*phi_k;
    phi_k_2_cos = cos(phi_k_2);
    phi_k_2_sin = sin(phi_k_2);
    delta_u_k = cu_sin*phi_k_2_sin + cu_cos*phi_k_2_cos; % Argument of lattitide correction
    delta_r_k = cr_sin*phi_k_2_sin + cr_cos*phi_k_2_cos; % Argument of radius correction
    delta_i_k = ci_sin*phi_k_2_sin + ci_cos*phi_k_2_cos; % Argument of inclination correction

    u_k = phi_k + delta_u_k; % Corrected argument of lattitude
    r_k = A*(1 - e*cos(E_k)) + delta_r_k; % Corrected argument of radius
    i_k = i_0 + delta_i_k + i_dot*t_k; % Corrected argument of inclination

    % Satellite position in orbital plane (x_k_dash, y_k_dash)
    x_k_dash = r_k * cos(u_k);
    y_k_dash = r_k * sin(u_k);

    % Corrected longitude of ascending node  % accounting for earth rotation rate
    omg_k = omg0 + (omg_dot - omg_e_dot)*t_k - omg_e_dot*t_0e;

    % Satellite Position in ECEF coordinates (x_k, y_k, z_k)
    x_k = x_k_dash*cos(omg_k) - y_k_dash*cos(i_k)*sin(omg_k);
    y_k = x_k_dash*sin(omg_k) + y_k_dash*cos(i_k)*cos(omg_k);
    z_k = y_k_dash*sin(i_k);
    
    P = [x_k; y_k; z_k]; % satellite position
    
end

