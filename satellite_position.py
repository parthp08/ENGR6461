import numpy as np

def sat_pos(eph_data, t_rcv, tau):
    """
    Calculate the Position of Satellite in ECEF frame given the ephemeris data of the satellite
    and GPS receiving time of the week

    Inputs
    -------
    eph_data : np.array, shape(21,), ephemeris data for a satellite
    t_rcv : float, receiving GPS time of the week, seconds
    tau: float, transmission time between sat and user, seconds

    Returns
    --------
    (x, y, z) : tuple, floats, Position of Satellte in ECEF frame, in meters
    
    References
    -----------
    [1] https://moodle.concordia.ca/moodle/pluginfile.php/3799910/mod_resource/content/1/Project_Appendix.pdf
    [2] https://moodle.concordia.ca/moodle/pluginfile.php/3799909/mod_resource/content/1/ENGR6461_project.pdf
    [3] 'Global Positioning System: Theory and Applications'; edited by Parkinson, 
         Spilker, Axelrad, Enge; AIAA; 1996
    [4] 'Aerospace Navigation Systems'; edited by A.V.Nebylv, J.Watson
    [5] https://ascelibrary.org/doi/pdf/10.1061/9780784411506.ap03
    [6] GNSS Applications and Methods; Scott Gleason, Demoz Gebre-Egzibher
    
    """

    ### eph parameters (Ref [2], page 2)
    cr_sin, delta_n, M_0, cu_cos, e, cu_sin, sqrta, t_0e, ci_cos, \
    omg0, ci_sin, i_0, cr_cos, w, omg_dot, i_dot, tgd, \
    toc, af2, af1, af0 = tuple(eph_data)

    ### Satellite Position Calculation Algorithm (Ref[1] page 3, Ref [3] page 138)
    mu = 3.986005e14    # m^3 / s^2 # Earth universal gravitational parameter
    omg_e_dot = 7.2921151467e-5 # Earth Rotation Rate
    a = sqrta**2    # semimajor axis
    
    n0 = np.sqrt(mu / a**3) # computed mean motion # rad/sec
    
    #tau = 0.075
    t = t_rcv - tau # t == t_tr

    t_k = t - t_0e  # time from ephemeris reference time(t_0e)
    # from the footnote(page 138 Ref [3])
    # accounting for beginning or end of week crossovers
    if t_k > 302400:  # seconds
        t_k = t_k - 604800
    if t_k < -302400:
        t_k = t_k + 604800
    
    n = n0 + delta_n    # corrected mean motion
    
    M_k = M_0 + n*t_k   # Mean anomaly
    
    pi = 3.1415926535898    # GPS standard value for pi
    
    # (Ref [1] page 3(Note))
    # E_k = M_k + e*np.sin(E_k)    # Kepler's Equ for the eccentric anomaly E_k # rad
    E_old = M_k
    E_new = M_k + e*np.sin(E_old)
    while E_new - E_old >= 1e-8:  # stop iteration when E_new - E_old is less than 1e-8
        E_old = E_new
        E_new = M_k + e*np.sin(E_old)
    E_k = E_new # Essentric anomaly at time t_k
    
    # True anomaly as a function of E_k at time t_k
    v_k = np.arctan2(np.sqrt(1 - e**2)*np.sin(E_k), np.cos(E_k) - e)

    E_k = np.arccos((e + np.cos(v_k)) / (1 + e*np.cos(v_k)))    # Eccentric anomaly

    phi_k = v_k + w   # Argument of latitude
    
    # Second harmonic perturbations
    phi_k_2 = 2*phi_k
    phi_k_2_cos = np.cos(phi_k_2)
    phi_k_2_sin = np.sin(phi_k_2)
    delta_u_k = cu_sin*phi_k_2_sin + cu_cos*phi_k_2_cos # Argument of lattitide correction
    delta_r_k = cr_sin*phi_k_2_sin + cr_cos*phi_k_2_cos # Argument of radius correction
    delta_i_k = ci_sin*phi_k_2_sin + ci_cos*phi_k_2_cos # Argument of inclination correction

    u_k = phi_k + delta_u_k # Corrected argument of lattitude
    r_k = a*(1 - e*np.cos(E_k)) + delta_i_k # Corrected argument of radius
    i_k = i_0 + delta_i_k + i_dot*t_k # Corrected argument of inclination

    # Satellite position in orbital plane (x_k_dash, y_k_dash)
    x_k_dash = r_k * np.cos(u_k)
    y_k_dash = r_k * np.sin(u_k)

    # Corrected longitude of ascending node  # accounting for earth rotation rate
    omg_k = omg0 + (omg_dot - omg_e_dot)*t_k - omg_e_dot*t_0e

    # Satellite Position in ECEF coordinates (x_k, y_k, z_k)
    x_k = x_k_dash*np.cos(omg_k) - y_k_dash*np.cos(i_k)*np.sin(omg_k)
    y_k = x_k_dash*np.sin(omg_k) + y_k_dash*np.cos(i_k)*np.cos(omg_k)
    z_k = y_k_dash*np.sin(i_k)

    # correction for earth rotation # Ref [6] page 61
    # by rotating position of sat(x_k,y_k,z_k) by small amount to obtain correct value of sat position
    x_k_n = np.cos(omg_e_dot*tau)*x_k + np.sin(omg_e_dot*tau)*y_k
    y_k_n = -np.sin(omg_e_dot*tau)*x_k + np.cos(omg_e_dot*tau)*y_k
    z_k_n = z_k

    return x_k_n, y_k_n, z_k_n


## Tests
if __name__ == "__main__":

    from scipy.io import loadmat
    data = loadmat('project_data.mat')
    # column i has information for satellite i
    # each row has 21 ephemeris info for satellites
    eph_data = data['eph']
    # print(eph_data.shape) # (21,6)

    eph_sat1 = eph_data[:, 1]
    # print(eph_sat1.shape)
    eph_sat1_tu = tuple(eph_sat1)
    # print(eph_sat1_tu)

    iono = data['iono'][0]
    t_rcv = iono[0]
    print(t_rcv)

    tau = 0.075 # initial estimation

    print(sat_pos(eph_sat1, t_rcv, tau))

