%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Least Square solution of GPS pseudorange equation
%
%%% references
% -------------
%[1] 'GNSS Applications and Methods', Scott Gleason, Demoz Gebre-Egzibher
%
%%% inputs
% ----------
% P_sat_arr : array, size(3,*), array of satellites' position in ECEF 
%             frame, meters
% Pu_0 : array, size(3,1), initial user position guess, meters
% pr_arr_data : array, size(1,*), psedorange data received from satellites,
%               meters
%
%%% outputs
%-----------
% delta_r : array, size(4,1),[delta_x, delta_y, delta_z, c*delta_b], meters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delta_r = least_square_sol(P_sat_arr, Pu_0, pr_arr_data)

    % number of satellite
    n = size(P_sat_arr, 2); % 2 for size of columns

    % vector dimension
    p = size(P_sat_arr(:,1), 1); % 1 for size of rows

    % initialize H matrix
    H = ones(n, p+1);
    
    % initialze rho_0
    rho_0 = ones(n,1);
    
    % fill rho_0 vector and H matrix
    for i = 1:n
        H(i,1:p) = (-unit_vector(Pu_0, P_sat_arr(:,i)));
        rho_0(i,1) = norm(P_sat_arr(:,i) - Pu_0);
    end
    
    H_T = H';
    H_terms = (inv(H_T*H))*(H_T);
    
    delta_rho = (pr_arr_data)' - rho_0;
    
    % delta_r = P_u - Pu_0
    delta_r = H_terms*delta_rho;
    
end