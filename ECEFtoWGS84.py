from math import sin, cos, atan, atan2, sqrt, pi, copysign
from numpy import rad2deg

"""
References
-----------

    [1] 'Aerospace Navigation Systems'; edited by A.V.Nebylv, J.Watson, page 32-33
"""

def ECEFtoWGS84(x,y,z):

    ## WGS84 Model
    a = 6378137.0  # meters # semi-major axis
    b = 6356752.314  # meters # semi-minor axis
    e = sqrt(1 - (b**2 / a**2))# eccentricity
    
    ### Step 1
    if (x**2 + y**2) == 0:  # if True algorithm stops here
        if z > 0:
            phi = pi/2
            h = z - b
        else: # z < 0
            phi = -pi/2
            h = -z - b
        return "Longitude cannot be defined", phi, h

    ### Step 2
    # assuming x**2 + y**2 > 0
    if x == 0:
        lambda_ = pi/2
    elif x > 0:
        lambda_ = atan(y/x)
    elif x < 0 and y >= 0 :
        lambda_ = pi + atan(y/x)
    elif x < 0 and y < 0:
        lambda_ = -pi + atan(y/x)
    
    ### Step 3
    # sensitivity levels
    e_phi = 0.01
    e_h = 0.01

    # starting values
    phi_s = pi/2
    d_s = a / sqrt(1 - e**2)
    h_s = 0

    ### Step 4
    # computing new values related to the latitude
    phi_n = atan2(z*(d_s+h_s) , ((sqrt(x**2 + y**2))*((d_s*(1 - e**2))+h_s)))
    d_n = a / sqrt(1 - ((e*sin(phi_n))**2))

    ### Step 5
    # computing new altitude value
    if phi_n <= pi/4 or phi_n >= -pi/4:
        h_n = (sqrt(x**2 + y**2)/cos(phi_n)) - d_n
    if phi_n > pi/4 or phi_n < -pi/4:
        h_n = (z/sin(phi_n)) - (d_n*(1 - e**2))

    ## Step 6
    # Evaluating the convergence of the obtained value for phi
    while abs(phi_n - phi_s) >= e_phi and abs(h_n - h_s) >= e_h: # if not converged
        phi_s = phi_n
        d_s = d_n
        h_s = h_n

        ### Step 4
        # computing new values related to the latitude
        phi_n = atan2(z*(d_s+h_s), (sqrt(x**2 + y**2))*(d_s*(1 - e**2) + h_s))
        d_n = a / sqrt(1 - ((e*sin(phi_n))**2))

        ### Step 5
        # computing new altitude value
        if phi_n <= pi/4 or phi_n >= -pi/4:
            h_n = (sqrt(x**2 + y**2)/cos(phi_n)) - d_n
        if phi_n > pi/4 or phi_n < -pi/4:
            h_n = (z/sin(phi_n)) - (d_n*(1 - e**2))

    # when converged
    # abs(phi_n - phi_s) >= e_phi and abs(h_n - h_s) >= e_h
    return rad2deg(lambda_), rad2deg(phi_n), h_n


## Tests
if __name__ == "__main__":

    print(ECEFtoWGS84(2.00, 2.00, 14412.42)) # 45, 89.997169, -6342339.81
    print(ECEFtoWGS84(15, 200, 13074.60)) # 85.710847, 89.794489, -6343677.35
    print(ECEFtoWGS84(0.00, 250.00, 10100.42)) # 90, 89.729935, -6346553.97
    # print(ECEFtoWGS84(1.00, 0.00, .00)) # 0.000000, 0.000000, -6378136.00
    print(ECEFtoWGS84(1.02, 0.00, 30446.61)) # 0.0,89.999203, -6326305.63

