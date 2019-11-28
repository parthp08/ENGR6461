%% plots of satellites position in ECEF frame
sat_pos_x = [P_a(1), P_b(1), P_c(1), P_d(1), P_e(1), P_f(1)];
sat_pos_y = [P_a(2), P_b(2), P_c(2), P_d(2), P_e(2), P_f(2)];
sat_pos_z = [P_a(3), P_b(3), P_c(3), P_d(3), P_e(3), P_f(3)];
plot3(sat_pos_x, sat_pos_y, sat_pos_z, 'o');

%% plots of user location with satellites position in ECEF frame
pos_x = [P_u(1), P_a(1), P_b(1), P_c(1), P_d(1), P_e(1), P_f(1)];
pos_y = [P_u(1), P_a(2), P_b(2), P_c(2), P_d(2), P_e(2), P_f(2)];
pos_z = [P_u(1), P_a(3), P_b(3), P_c(3), P_d(3), P_e(3), P_f(3)];
plot3(pos_x, pos_y, pos_z, 'o');

%% On 3D Earth
sat_pos_x = [P_a(1), P_b(1), P_c(1), P_d(1), P_e(1), P_f(1)];
sat_pos_y = [P_a(2), P_b(2), P_c(2), P_d(2), P_e(2), P_f(2)];
sat_pos_z = [P_a(3), P_b(3), P_c(3), P_d(3), P_e(3), P_f(3)];
pos_x = [P_u(1), P_a(1), P_b(1), P_c(1), P_d(1), P_e(1), P_f(1)];
pos_y = [P_u(1), P_a(2), P_b(2), P_c(2), P_d(2), P_e(2), P_f(2)];
pos_z = [P_u(1), P_a(3), P_b(3), P_c(3), P_d(3), P_e(3), P_f(3)];
earth_sphere('m')
hold on
plot3(pos_x, pos_y, pos_z, 'o')
%plot3(sat_pos_x, sat_pos_y, sat_pos_z, 'o');
