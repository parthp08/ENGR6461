import numpy as np

## TODO: Write Docstring

def iono_error(iono_data, lat_u, lon_u, E_s, A_s):
    """
    lat_u = users apporx geodetic latitude
    lon_u = users apporx geodetic longitude
    E_s = elevation angle to (each) GPS Satellite
    A_s = Azimuth angle to (each) GPS Satellite

    """

    phi_u = lat_u
    lamda_u = lon_u

    # Ref [2] page 2
    tow, a0, a1, a2, a3, b0, b1, b2, b3 = tuple(iono_data)

    ### Klobuchar model for ionosphere  # Ref [1] page 4-6 and Ref [2] page 146-147 and Ref [4] page 78-79
    # Earth-centered angle between user and IP(ionosphere point # Ref[4] page 78)
    psi_ip = (0.0137/(0.11 + E_s)) - 0.022 # semicircles
    
    # subionospheric latitide
    phi_ip = phi_u + psi_ip*np.cos(A_s)    # semicircles
    if phi_ip > 0.416:
        phi_ip = 0.416
    if phi_ip < -0.416:
        phi_ip = -0.416
    
    # subionospheric longitude
    lamda_ip = lamda_u + psi_ip*np.sin(A_s)/np.cos(phi_ip)  # semicircles

    # geomagnetic lattitude of the subionospheric location looking toward each GPS satellite
    phi_m = phi_ip + 0.064*np.cos(lamda_ip-1.617)   # semicircles

    # local time(t) at the subionospheric point
    t_gps = tow     # is it same ????????????????????????????????????????
    t_ip = 4.32e4*lamda_ip + t_gps # seconds       # t_gps == ???????????????????????
    if t_ip > 86400:
        t_ip = t_ip - 86400
    if t_ip < 0:
        t_ip = t_ip + 86400
    
    # obliquity factor / slant factor (dimensionless)
    # to convert to slant time delay 
    F_ip = 1.0 + 16.0*((0.53-E_s)**3)

    # Amplitude # in seconds
    AM_ip = 0
    a = [a0,a1,a2,a3]
    for i in range(len(a)):
        AM_ip += a[i]*(phi_m**i)
    if AM_ip < 0:
        AM_ip = 0

    # Period # in seconds
    PER_ip = 0
    b = [b0,b1,b2,b3]
    for i in range(len(b)):
        PER_ip += b[i]*(phi_m**i)
    if PER_ip < 72000:
        PER_ip = 72000

    # Phase # in rad
    X_ip = 2*np.pi*(t_ip - 50400) / PER_ip

    # Ionospheric Time delay (T_iono)  # in seconds
    if np.abs(X_ip) >= 1.57:
        T_iono = F_ip*5e-9
    else:   # for np.abs(X_ip) < 1.57
        T_iono = F_ip*(5e-9 + AM_ip*(1 - (X_ip**2)/2 + (X_ip**4)/24))

    # Ionosphere Error # in meters
    c = 299792458   # speed of light(GPS constant) # m/s
    I_e = c * T_iono    # m
    
    return I_e



## Tests
if __name__ == "__main__":

    from scipy.io import loadmat
    data = loadmat('project_data.mat')
    """
        References
        -----------
        [1] https://moodle.concordia.ca/moodle/pluginfile.php/3799910/mod_resource/content/1/Project_Appendix.pdf
        [2] https://moodle.concordia.ca/moodle/pluginfile.php/3799909/mod_resource/content/1/ENGR6461_project.pdf
        [3] 'Global Positioning System: Theory and Applications'; edited by Parkinson, 
            Spilker, Axelrad, Enge; AIAA; 1996
        [4] 'Aerospace Navigation Systems'; edited by A.V.Nebylv, J.Watson

    """
    iono = data['iono'].reshape(9,)

    print(iono_error(iono, -15, 15, 10, 20))  # Just to check if this works # might not be valid inputs 
    ## IS this can be -ve ??

