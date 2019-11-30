%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute the tropospheric errors in pseudorange mesurement
%
%%% references
% -------------
%[1] https://moodle.concordia.ca/moodle/pluginfile.php/3799910/mod_resource/content/1/Project_Appendix.pdf
%[2] 'Aerospace Navigation Systems'; edited by A.V.Nebylv, J.Watson
%
%%% inputs
% ----------
% T_amb : float, ambient air tempreature, deg celsius 
% P_amb : float, ambient air pressure, kPa 
% P_vap : float, ambient vapour pressure, kPa
% E_s : float, elevation of sat from the user local tangent plane, rad
%
%%% outputs
%-----------
% T_e : float, tropospheric delay, meters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T_e = tropo_error(T_amb, P_amb, P_vap, E_s)

    % zenith delay of the dry component
    Kd = (1.55208e-4*P_amb*(40136.0 + 148.72*T_amb)) / (T_amb+273.16);   % m
    
    % zenith delay of the wet component
    Kw = (-0.282*P_vap/(T_amb+273.16)) + (8307.2*P_vap/((T_amb+273.16)^2)); % m

    % Troposhperic delay correction (T_e), m
    temp1 = sqrt(E_s*E_s + 1.904e-3);
    temp2 = sqrt(E_s*E_s + 0.6854e-3);
    T_e = (Kd/sin(temp1)) + (Kw/sin(temp2));   % m

end

