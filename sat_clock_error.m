%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute satellite clock correction
%
%%% references
% -------------
%[1] https://moodle.concordia.ca/moodle/pluginfile.php/3799910/mod_resource/content/1/Project_Appendix.pdf
%[2] https://moodle.concordia.ca/moodle/pluginfile.php/3799909/mod_resource/content/1/ENGR6461_project.pdf
%[3] 'Global Positioning System: Theory and Applications'; edited by 
%     Parkinson, Spilker, Axelrad, Enge; AIAA; 1996
%[4] 'Aerospace Navigation Systems'; edited by A.V.Nebylv, J.Watson
%[5] https://ascelibrary.org/doi/pdf/10.1061/9780784411506.ap03
%
%%% inputs
% ----------
% eph_data : array, size(21,*), ephemeris data received from the GPS
%            satellites
% t_rcv : flot, float, receiving GPS time of the week, seconds
% tau : float, transmission time between sat and user, seconds
%
%%% outputs
%-----------
% delta_t_sv : array, size(1,*), float, satellite clock correction, seconds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delta_t_sv = sat_clock_error(eph_data, t_rcv, tau)

    % eph parameters (Ref [2], page 2)
    %cr_sin = eph_data(1,:);
    delta_n = eph_data(2,:) * pi;
    M_0 = eph_data(3,:) * pi;
    %cu_cos = eph_data(4,:);
    e = eph_data(5,:);
    %cu_sin = eph_data(6,:);
    sqrta = eph_data(7,:);
    t_0e = eph_data(8,:);
    %ci_cos = eph_data(9,:);
    %omg0 = eph_data(10,:) * pi;
    %ci_sin = eph_data(11,:);
    %i_0 = eph_data(12,:) * pi;
    %cr_cos = eph_data(13,:);
    %w = eph_data(14,:) * pi;
    %omg_dot = eph_data(15,:) * pi;
    %i_dot = eph_data(16,:) * pi;
    tgd = eph_data(17,:);
    toc = eph_data(18,:);
    af2 = eph_data(19,:);
    af1 = eph_data(20,:);
    af0 = eph_data(21,:);

    F = -4.442807633e-10;    % relativistic constant

    % inividual satellite time corrected to GPS system time 't'
    % trasmission time  % Ref [2] page 1
    % tau = 0.075 % initial estimation    % tau = t_rcv - t_tr
    t_tr = t_rcv - tau; % t == t_tr

    t_k = t_tr - t_0e;  % time from ephemeris reference time(t_0e)
    % account for end of the week crossover correction % Ref [1] page 2
    for i = 1:length(t_k)
        if t_k(i) > 302400  % seconds
            t_k(i) = t_k(i) - 604800;
        end
        if t_k(i) < -302400
            t_k(i) = t_k(i) + 604800;
        end
    end
        
    mu = 3.986005e14;    % m^3 / s^2 % Earth universal gravitational parameter
    a = sqrta.^2;    % semimajor axis

    n0 = sqrt(mu ./ a.^3); % computed mean motion % rad/sec

    n = n0 + delta_n;    % corrected mean motion
    
    M_k = M_0 + n.*t_k;   % Mean anomaly
    
    % (Ref [1] page 3(Note))
    % E_k = M_k + e*np.sin(E_k)    % Kepler's Equ for the eccentric anomaly E_k % rad
    E_old = M_k;
    E_new = M_k + e.*sin(E_old);
    for i = 1:length(E_new)
        while E_new(i) - E_old(i) >= 1e-8  % stop iteration when E_new - E_old is less than 1e-8
            E_old(i) = E_new(i);
            E_new(i) = M_k(i) + e(i).*sin(E_old(i));
        end
    end
    E_k = E_new; % Essentric anomaly at time t_k
    
    % Relativistic correction term
    T_rel = F.*e.*sqrta.*sin(E_k);   % sec

    % Ref [1] page 1 and Ref [5]
    delta_t_sv = af0 + af1.*(t_tr - toc) + af2.*((t_tr - toc).^2) + T_rel - tgd;
    % t = t_tr - delta_t_sv
    
end

