r""" This script is the inner workings of the Propagate Project

This code will simplify the functions for the user
It will connect the functions into one simple code
Allows for Text File reading
Includes final conversions to degrees
Prints Results to Text File

Author: Thomas J Susi
"""

import numpy as np
from astro import RV2COE
from astro import PROPAGATE
from astro import constants
import pdb

filename = input("Name File to Interpret: ")
with open(filename) as file, open("PROPAGATE_Results.txt", "w") as finalfile:

    print ("For Inputting Lines: Line n = n-1")
    start = input("First Line to Evaluate: "),
    end = input("Last Line to Evaluate: "),
    mu = float(input("Input mu Value: "))

    for count in range(int(start[0]), int(end[0])):

        line = open (filename, "r").readlines()[count].split()
        Line_0 = ['Line {}\n\t{}\n'.format(count+1,line)]

        r_i = float(line[0])
        r_j = float(line[1])
        r_k = float(line[2])
        v_i = float(line[3])
        v_j = float(line[4])
        v_k = float(line[5])

        change_t = float(line[6]) * 60

        r_input = np.array([r_i, r_j, r_k])
        v_input = np.array([v_i, v_j, v_k])

        Line_1 = ("The mu Constant: \n\t", str(mu), "\n", 'Original Radius {} (km)\n\t{}\n'.format(count+1,r_input))
        Line_2 = ('Original Velocity {} (km/s)\n\t{}\n'.format(count+1,v_input))
    
        a_i, e_i, e_vector, i_i, raan_i, w_i, theta_i, p, r_p, r_a, n_vector, period = RV2COE.RV2COE(r_input, v_input, mu)

        a_f, e_f, i_f, raan_f, w_f, theta_f = PROPAGATE.update(a_i, e_i, i_i, raan_i, w_i, theta_i, change_t, mu)

        r_eci, v_eci = PROPAGATE.COE2RV(a_f, e_f, i_f, raan_f, w_f, theta_f, mu)

        r"""Print The Results in Full"""
        Line_3 = ["Semi Major Axis in Kilometers:\n\t", str(a_f), "\n"]
        Line_4 = ["Eccentricity:\n\t", str(e_f), "\n"]
        Line_5 = ["Eccentricity Vector:\n\t", str(e_vector), "\n"]
        Line_6 = ["Initial Inclination in Degrees:\n\t", str(np.rad2deg(i_i)), "\n"]
        Line_7 = ["Initial RAAN in Degrees:\n\t", str(np.rad2deg(raan_i)), "\n"]
        Line_8 = ["Initial Argument of Periapsis in Degrees:\n\t", str(np.rad2deg(w_i)), "\n"]
        Line_9 = ["Initial True Anomaly in Degrees:\n\t", str(np.rad2deg(theta_i)), "\n"]
        Line_61 = ["Final Inclination in Degrees:\n\t", str(np.rad2deg(i_f)), "\n"]
        Line_71 = ["Final RAAN in Degrees:\n\t", str(np.rad2deg(raan_f)), "\n"]
        Line_81 = ["Final Argument of Periapsis in Degrees:\n\t", str(np.rad2deg(w_f)), "\n"]
        Line_91 = ["Final True Anomaly in Degrees:\n\t", str(np.rad2deg(theta_f)), "\n"]
        Line_10 = ["Semi Latus Rectum in Kilometers:\n\t", str(p), "\n"]
        Line_11 = ["Radius of Periapsis in Kilometers:\n\t", str(r_p), "\n"]
        Line_12 = ["Radius of Apoapsis in Kilometers:\n\t", str(r_a), "\n"]
        Line_14 = ["Line of Nodes Vector in Kilometers:\n\t", str(n_vector), "\n"]
        Line_15 = ["Orbit Period in Seconds:\n\t", str(period), "\n"]
        Line_16 = ["New Radius Vector (ECI):\n", str(r_eci), "\n"]
        Line_17 = ["New Velocity Vector (ECI):\n", str(v_eci), "\n"]
        Line_18 = ["\n\n\n"]

        finalfile.writelines(Line_0)
        finalfile.writelines(Line_1)
        finalfile.writelines(Line_2)
        finalfile.writelines(Line_3)
        finalfile.writelines(Line_4)
        finalfile.writelines(Line_5)
        finalfile.writelines(Line_6)
        finalfile.writelines(Line_7)
        finalfile.writelines(Line_8)
        finalfile.writelines(Line_9)
        finalfile.writelines(Line_61)
        finalfile.writelines(Line_71)
        finalfile.writelines(Line_81)
        finalfile.writelines(Line_91)
        finalfile.writelines(Line_10)
        finalfile.writelines(Line_11)
        finalfile.writelines(Line_12)
        finalfile.writelines(Line_14)
        finalfile.writelines(Line_15)
        finalfile.writelines(Line_16)
        finalfile.writelines(Line_17)
        finalfile.writelines(Line_18)
finalfile.close()

print ("Complete")
