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
from LonLatToElvAz import LonLatToElvAz

# constants 
c = 299792458   # m/sec

# initial user position and bias
xu0, yu0, zu0 = (0, 0, 0)
b = 0

# receiving GPS time of the week
t_rcv = iono[0]
tau = 0.075 # intial estimation of transmission time between sat and user

# ionospheric error
I_e = 0 # initially set to zero

# tropospheric error
T_e = 0 # initally set to zero

# initial guess
delta_x, delta_y, delta_z, delta_b = 100, 100, 100, 100

while norm(delta_x, delta_y, delta_z) > 1e-3 and delta_b > 1:

    delta_rho_arr = [] # delta_rho for 6 satellite
    sat_pos_x = []  # positions of all 6 sat in x-direction
    sat_pos_y = []
    sat_pos_z = []
    e_array = []    # unit vectors coordinates for 6 satellites

    for i in range(0,6):    # 6 satellites
        # estimate satellite clock error
        delta_t_sv = sat_clock_error(eph_data[:,i], t_rcv, tau)

        # get satellite position
        xs, ys, zs = sat_pos(eph_data[:,i], t_rcv, tau)
        sat_pos_x.append(xs)
        sat_pos_y.append(ys)
        sat_pos_z.append(zs)

        # pseudo-range for sat
        rho = pr_data[i]

 #       -----
        ## ?? where do i calulate these so it can be added to pseudorange corrections
        # data for calculating iono_error and tropo_error
        T_amb = 15  # deg C
        P_amb = 101.325 # kPa
        P_vap = 0.85    # kPa
        lon_u, lat_u, h_u = ECEFtoWGS84(xu0,yu0,zu0)
        # E should be in degree
        E_s, A_s = LonLatToElvAz(xu0, yu0, zu0, xs, ys, zs)

        # Ionosphere Corrections
        I_e = iono_error(iono, lat_u, lon_u, E_s, A_s)

        # Troposhpere Corretions
        T_e = tropo_error(T_amb, P_amb, P_vap, E_s)
#        -----

        # apply correction to pseudo-range
        # new_rho =  pseudorange - sat_clock_error - user_clock_error - iono_error - tropo_error
        rho_corrected = rho - c*delta_t_sv - c*b - I_e - T_e  #?????# not sure of this equation

        # new signal trasmission time (tau)
        tau = rho_corrected / c

        # while loop for convergence 
        # unit vector from user to sat # now from origin
        e = unit_vect_u2s(xu0,yu0,zu0,xs,ys,zs) 
        e_array.append(e)

        rho_0 = vec_magnitude(xu0,yu0,zu0,xs,ys,zs)
        rho = rho_corrected    # pseudo-range for sat # corrected
        delta_rho_arr.append(rho - rho_0)

    H = H_matrix(e_array)
    H_terms = H_terms_sol(H)

    # Least Square solution
    delta_r = np.matmul(H_terms, np.array(delta_rho_arr).reshape(6,1))
    delta_x, delta_y, delta_z, delta_b = delta_r.reshape(4,)

    # pos_new = pos_old + delta_pos # GNSS book page 61 # Aero Nav book page 36
    xu0 = xu0 + delta_x
    yu0 = yu0 + delta_y
    zu0 = zu0 + delta_z
    b = b + delta_b

    # ## ?? where do i calulate these so it can be added to pseudorange corrections
    # # data for calculating iono_error and tropo_error
    # T_amb = 15  # deg C
    # P_amb = 101.325 # kPa
    # P_vap = 0.85    # kPa
    # lon_u, lat_u, h_u = ECEFtoWGS84(xu0,yu0,zu0)
    # # E should be in degree
    # E_s, A_s = LonLatToElvAz(xu0, yu0, zu0, xs, ys, zs)

    # # Ionosphere Corrections
    # I_e = iono_error(iono, lat_u, lon_u, E_s, A_s)

    # # Troposhpere Corretions
    # T_e = tropo_error(T_amb, P_amb, P_vap, E_s)


# after solution converges
xu = xu0
yu = yu0
zu = zu0
bu = b
print(xu, yu, zu, bu)


# # plot the users and satellites position on the map
# # Ref : https://matplotlib.org/2.1.1/gallery/mplot3d/scatter3d.html
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# # points to be ploteed 
# # user location and sat positions
# x_arr = [xu] + sat_pos_x
# y_arr = [yu] + sat_pos_y
# z_arr = [zu] + sat_pos_z

# ax.scatter(sat_pos_x, sat_pos_y, sat_pos_z) # sat positions
# ax.scatter(xu, yu, zu)   # user positions
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')

# plt.show()
### Plotting is not working after addition of Iono and Tropo Errors ?????


