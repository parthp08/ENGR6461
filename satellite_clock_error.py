import numpy as np

"""
    References
    -----------
    [1] https://moodle.concordia.ca/moodle/pluginfile.php/3799910/mod_resource/content/1/Project_Appendix.pdf
    [2] https://moodle.concordia.ca/moodle/pluginfile.php/3799909/mod_resource/content/1/ENGR6461_project.pdf
    [3] 'Global Positioning System: Theory and Applications'; edited by Parkinson, 
         Spilker, Axelrad, Enge; AIAA; 1996
    [4] 'Aerospace Navigation Systems'; edited by A.V.Nebylv, J.Watson
    [5] https://ascelibrary.org/doi/pdf/10.1061/9780784411506.ap03

"""

def sat_clock_error(eph_data, t_rcv, tau):
    """
    Estimate satellite clock bias

    # Ref [1] page 2
    # Ref [5]
     
    """

    ### eph parameters (Ref [2], page 2)
    cr_sin, delta_n, M_0, cu_cos, e, cu_sin, sqrta, t_0e, ci_cos, \
    omg0, ci_sin, i_0, cr_cos, w, omg_dot, i_dot, tgd, \
    toc, af2, af1, af0 = tuple(eph_data)

    F = -4.442807633e-10    # relativistic constant

    # inividual satellite time corrected to GPS system time 't'
    # trasmission time  # Ref [2] page 1
    # tau = 0.075 # initial estimation    # tau = t_rcv - t_tr
    t_tr = t_rcv - tau # t == t_tr

    t_k = t_tr - t_0e  # time from ephemeris reference time(t_0e)
    # account for end of the week crossover correction # Ref [1] page 2
    if t_k > 302400:  # seconds
        t_k = t_k - 604800
    if t_k < -302400:
        t_k = t_k + 604800

    mu = 3.986005e14    # m^3 / s^2 # Earth universal gravitational parameter
    a = sqrta**2    # semimajor axis

    n0 = np.sqrt(mu / a**3) # computed mean motion # rad/sec

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

    # Relativistic correction term
    T_rel = F*e*sqrta*np.sin(E_k)   # sec

    # Ref [1] page 1 and Ref [5]
    delta_t_sv = af0 + af1*(t_tr - toc) + af2*((t_tr - toc)**2) + T_rel - tgd
    # t = t_tr - delta_t_sv

    return delta_t_sv


# Tests
if __name__ == "__main__":
    from scipy.io import loadmat
    data = loadmat('project_data.mat')
    eph_data = data['eph']
    eph1 = eph_data[:,1]
    t_rcv = data['iono'][0][0]
    print(sat_clock_error(eph1, t_rcv))

