"""
Ref : Aerospace Navigation, page 36-37
"""

import numpy as np

from scipy.io import loadmat
data = loadmat('project_data.mat')
eph_data = data['eph']
iono = data['iono']
iono = iono.reshape(9,)
pr_data = data['pr']
pr_data = pr_data.reshape(6,)

from ECEFtoWGS84 import ECEFtoWGS84
from satellite_position import sat_pos
from satellite_clock_error import sat_clock_error
from ionospheric_error import iono_error
from tropospheric_error import tropo_error
from unit_vectors_u2s import unit_vect_u2s
from vector_magnitude import vec_magnitude
from norm import norm
from H_matrix import H_matrix, H_terms_sol

def step2to4(xu, yu, zu, delta_t_u, t_rcv, tau, c):

    delta_rho_arr = [] # delta_rho for 6 satellite
    sat_pos_x = []  # positions of all 6 sat in x-direction
    sat_pos_y = []
    sat_pos_z = []
    e_array = []    # unit vectors coordinates for 6 satellites

    for i in range(0,6):
        
        ### Step: 2
        Xs, Ys, Zs = sat_pos(eph_data[:,i], t_rcv, tau)
        sat_pos_x.append(Xs)
        sat_pos_y.append(Ys)
        sat_pos_z.append(Zs)
        r_i = vec_magnitude(Xs, Ys, Zs, xu, yu, zu)
        rho_ = r_i #+ c*delta_t_u

        e = unit_vect_u2s(xu,yu,zu,Xs,Ys,Zs) 
        e_array.append(e)

        delta_rho_arr.append(rho_ - pr_data[i])


    H = H_matrix(e_array)
    H_terms = H_terms_sol(H)

    # Least Square solution
    delta_r = np.matmul(H_terms, np.array(delta_rho_arr).reshape(6,1))
    delta_x, delta_y, delta_z, delta_t_u_ = delta_r.reshape(4,)
    delta_cb_u = -c*delta_t_u_
    # print(delta_x, delta_y, delta_z, delta_cb_u)

    xu = xu + delta_x
    yu = yu + delta_y
    zu = zu + delta_z
    delta_t_u = delta_t_u + delta_t_u_

    return xu, yu, zu, delta_x, delta_y, delta_z, delta_t_u


# constants 
c = 299792458   # m/sec

# receiving GPS time of the week
t_rcv = iono[0]
tau = 0.075 # intial estimation of transmission time between sat and user

### Step: 1
e_ = 1e-3   # sensitivity level
# starting values
xu, yu, zu, delta_t_u = 0, 0, 0, 0

### Step 2 to 4
xu, yu, zu, delta_x, delta_y, delta_z, delta_t_u = step2to4(xu, yu, zu, delta_t_u, t_rcv, tau, c)

### Step 5
while np.sqrt(delta_x**2 + delta_y**2 + delta_z**2) > e_:
    xu, yu, zu, delta_x, delta_y, delta_z, delta_t_u = step2to4(xu, yu, zu, delta_t_u, t_rcv, tau, c)

print(xu, yu, zu, delta_t_u)
