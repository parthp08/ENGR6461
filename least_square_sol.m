function [P_u, cb_u, delta_r] = least_square_sol(P_sat_arr, Pu_0, pr_arr_data)

    % number of satellite
    n = size(P_sat_arr, 2); % 2 for size of columns

    % vector dimension
    p = size(P_sat_arr(:,1), 1); % 1 for size of rows

    % initialize H matrix
    H = ones(n, p+1);
    
    % initialze rho_0
    rho_0 = ones(n,1);
    
    for i = 1:n
        H(i,1:p) = (-unit_vector(Pu_0, P_sat_arr(:,i)));
        rho_0(i,1) = norm(P_sat_arr(:,i) - Pu_0);
    end
    
    H_T = H';
    H_terms = (inv(H_T*H))*(H_T);
    
    delta_rho = (pr_arr_data)' - rho_0;
    
    % delta_r = P_u - Pu_0
    delta_r = H_terms*delta_rho;
    
    P_u = delta_r(1:3) + Pu_0;
    cb_u = delta_r(end);
    
end