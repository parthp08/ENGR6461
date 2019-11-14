import numpy as np

def unit_vect_u2s(xu,yu,zu,xs,ys,zs):
    """
    user position
    sat position

    returns components of the unit vector joining two points

    """

    # coordinates of Vector between two points
    V = [(xs-xu), (ys-yu), (zs-zu)]

    # distance between two points
    d = np.sqrt((xs-xu)**2 + (ys-yu)**2 + (zs-zu)**2)

    # coordinates of unit vector between two points
    e = [V[0]/d, V[1]/d, V[2]/d]

    return e
