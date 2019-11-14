## main calculations

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
from H_matrix import H_matrix, H_terms

#?? how to find user position
#?? then how to use errors in the solutions


# Trying Least Square Solution 😅😅
## delta_r = (H^T H)-1 H^T delta_rho
# assuming initial conditions to zero xu0,yu0,zu0 = (0,0,0) # origin
# rho_0 = from origin to sat_pos
# rho = from data
# delta_rho = rho - rho_0
# delta_r = r - r_o  == r   # r_0 == origin
# r => [xu,yu,zu,bcu]
# H-matrix elements => unit vectors # sat_pos from origin

xu0, yu0, zu0 = (0.0, 0.0, 0.0) # initial position # origin
t_rcv = iono[0] # receiving GPS time of the week

delta_rho_arr = [] # delta_rho for 6 satellite
sat_pos_x = []  # positions of all 6 sat in x-direction
sat_pos_y = []
sat_pos_z = []
e_array = []    # unit vectors coordinates for 6 satellites
for i in range(0,6):    # 6 satellites
    xs, ys, zs = sat_pos(eph_data[:,i], t_rcv) # sat position
    e = unit_vect_u2s(xu0,yu0,zu0,xs,ys,zs) # unit vector from user to sat # now from origin
    rho_0 = vec_magnitude(xu0,yu0,zu0,xs,ys,zs)
    rho = pr_data[i]    # pseudo-range for sat
    delta_rho_arr.append(rho - rho_0)
    e_array.append(e)
    sat_pos_x.append(xs)
    sat_pos_y.append(ys)
    sat_pos_z.append(zs)


H = H_matrix(e_array)

H_terms = H_terms(H)

# Least Square solution
delta_r = np.matmul(H_terms, np.array(delta_rho_arr).reshape(6,1))
print(delta_r)
xu, yu, zu, _ = delta_r.reshape(4,)

# plot the users and satellites position on the map
# Ref : https://matplotlib.org/2.1.1/gallery/mplot3d/scatter3d.html
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# points to be ploteed 
# user location and sat positions
x_arr = [xu] + sat_pos_x
y_arr = [yu] + sat_pos_y
z_arr = [zu] + sat_pos_z

ax.scatter(sat_pos_x, sat_pos_y, sat_pos_z) # sat positions
ax.scatter(xu, yu, zu, c='r')   # user positions
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()

