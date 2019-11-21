import numpy as np
from ECEFtoWGS84 import ECEFtoWGS84
from norm import norm

def LonLatToElvAz(xu,yu,zu,xs,ys,zs):

    lon, lat, h = ECEFtoWGS84(xu,yu,zu)
    X = np.array([xs-xu, ys-yu, zs-zu]).reshape(3,1)
    enu = Rotxyz2enu(X, lat, lon)

    A_s = np.arctan2(enu[0], enu[1])
    E_s = np.arcsin(enu[2] / norm(enu[0], enu[1], enu[2]))

    return E_s, A_s

    # return Elv, Az
    
def Rotxyz2enu(X, lat, lon):
    """
    https://github.com/kwogger/Autonomous-Mobile-Robotics/blob/master/matlab/ME597-5-GDOP/rotxyz2enu.m
    """

    R1 = rot((np.pi/2)+lon, 3)
    R2 = rot((np.pi/2)-lat, 1)

    R = np.matmul(R2, R1)

    enu = np.matmul(R, X)

    return enu


def rot(angle, axis):
    """
    https://github.com/kwogger/Autonomous-Mobile-Robotics/blob/master/matlab/ME597-5-GDOP/rot.m
    """

    R = np.eye(3)
    cang = np.cos(angle)
    sang = np.sin(angle)

    if axis == 1:
        R[1,1] = cang
        R[2,2] = cang
        R[1,2] = sang
        R[2,1] = -sang

    if axis == 2:
        R[0,0] = cang
        R[2,2] = cang
        R[0,2] = -sang
        R[2,0] = sang
        
    if axis == 3:
        R[0,0] = cang
        R[1,1] = cang
        R[1,0] = -sang
        R[0,1] = sang
        
    return R

