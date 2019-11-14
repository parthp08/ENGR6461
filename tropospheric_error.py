import numpy as np

## TODO: Write Docstring

def tropo_error(T_amb, P_amb, P_vap, E_s):
    """
    Hopfield's model # Ref [1] page 4 and Ref[4] page 81

    # Ref[1] page 4
    values and units of T_amb, P_amb, P_vap

    E_s = satellite Elevation angle  # in rad only

    """

    # zenith delay of the dry component
    Kd = (1.55208e-4*P_amb*(40136.0 + 148.72*T_amb)) / (T_amb+273.16)   # m
    
    # zenith delay of the wet component
    Kw = (-0.282*P_vap/(T_amb+273.16)) + (8307.2*P_vap/((T_amb+273.16)**2)) # m

    # Troposhperic delay correction (T_e), m
    temp1 = np.sqrt(E_s*E_s + 1.904e-3)
    temp2 = np.sqrt(E_s*E_s + 0.6854e-3)
    T_e = (Kd/np.sin(temp1)) + (Kw/np.sin(temp2))   # m

    return T_e


## Test
if __name__ == "__main__":

    # Ref [1] page 4 # using standard atmosphere values and Elevation of 12.86 deg
    T_amb = 15  # deg C
    P_amb = 101.325 # kPa
    P_vap = 0.85    # kPa
    E_s = np.deg2rad(12.86) # rad

    print(tropo_error(T_amb,P_amb,P_vap,E_s)) # should be 10.575 m

    """
        References
        -----------
        [1] https://moodle.concordia.ca/moodle/pluginfile.php/3799910/mod_resource/content/1/Project_Appendix.pdf
        [2] https://moodle.concordia.ca/moodle/pluginfile.php/3799909/mod_resource/content/1/ENGR6461_project.pdf
        [3] 'Global Positioning System: Theory and Applications'; edited by Parkinson, 
            Spilker, Axelrad, Enge; AIAA; 1996
        [4] 'Aerospace Navigation Systems'; edited by A.V.Nebylv, J.Watson

    """
