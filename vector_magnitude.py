import numpy as np

def vec_magnitude(x1,y1,z1,x2,y2,z2):
    """
    find the distance between two points (magnitude of the vector joining them)
    (Mangitude is always positive)
    """

    D = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

    return D

