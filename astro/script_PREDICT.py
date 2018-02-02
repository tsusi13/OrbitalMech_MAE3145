r""" This script is the inner workings of the PREDICT Project

This code will simplify the functions for the user
It will connect the functions into one simple code
Allows for Text File reading
Includes final conversions to degrees
Prints Results to Text File

Author: Thomas J Susi
"""

import numpy as np
from astro import PREDICT
from astro import COMFIX
from astro import PROPAGATE
from astro import RV2COE
from astro import time
from astro import planets
from astro import tle
import pdb

filename = input("Name File to Interpret: ")
with open(filename) as file, open("PREDICT_Results.txt", "w") as finalfile:

    print ("Input Observation Site Information")
    lat = np.deg2rad(float(input("Input Latitude(deg): ")))
    lon = np.deg2rad(float(input("Input Longitude(deg): ")))
    alt = float(input("Input Altitude(km): "))

    r_site_ecef = COMFIX.lla2ecef(lat, lon, alt)

    print ("Input Observation Window Information:")
    print ("Initial Time:")
    year_i = int(input("Input Initial Year: "))
    month_i = int(input("Input Initial Month: "))
    day_i = int(input("Input Initial Day: "))
    hour_i = int(input("Input Initial Hour: "))
    min_i = int(input("Input Initial Minute: "))
    sec_i = int(input("Input Initial Second: "))

    JD_i = time.date2jd(year_i, month_i, day_i, hour_i, min_i, sec_i)[0]

    print ("Final Time:")
    year_f = int(input("Input Final Year: "))
    month_f = int(input("Input Final Month: "))
    day_f = int(input("Input Final Day: "))
    hour_f = int(input("Input Final Hour: "))
    min_f = int(input("Input Final Minute: "))
    sec_f = int(input("Input Final Second: "))

    JD_f = time.date2jd(year_f, month_f, day_f, hour_f, min_f, sec_f)[0]

    time_stepmins = float(input("Input Time Step(mins): "))
    time_stepJD = (time_stepmins * 60) / 86400
    time_stepsec = time_stepmins * 60
    mu = 398600.5

    Line_ob = ['Observation Site Information:\n\tLatitude(deg): {}\n\tLongitude(deg): {}\n\tAltitude(km): {}\n'.format(np.rad2deg(lat), np.rad2deg(lon), alt)]
    finalfile.writelines(Line_ob)
    Line_times = ['Observation Time:\n\tJD Initial: {}\n\tJD Final: {}\n'.format(JD_i, JD_f)]
    finalfile.writelines(Line_times)

    print ("Choose TLES:")
    start_tle = int(input("First TLE(0 being first TLE): "))
    end_tle = int(input("Last TLE: "))
    century = int(input("Century of TLEs: "))

    tles = tle.get_tle(filename)

    for count in range(start_tle, end_tle):

        Line_0 = ['TLE {} for {}:\n'.format(count, tles[count].satname)]
        finalfile.writelines(Line_0)

        """TLE Initial Info Sorted"""
        i_0 = np.deg2rad(tles[count].inc)
        raan_0 = np.deg2rad(tles[count].raan)
        e_0 = tles[count].ecc
        w_0 = np.deg2rad(tles[count].argp)
        M_0 = np.deg2rad(tles[count].ma)
        n_0 = tles[count].mean_motion * 2 * np.pi / 86400
        nrate_0 = tles[count].ndot_over_2 * 2 * np.pi / 7464960000

        a_0 = (mu / (n_0**2))**(1/3)
        p_0 = a_0 * (1 - (e_0**2))

        epoch_year = tles[count].epoch_year + century
        epoch_dayfrac = tles[count].epoch_day

        epoch_mon, epoch_day, epoch_hr, epoch_min, epoch_sec = time.dayofyr2mdhms(epoch_year, epoch_dayfrac)
        JD_epoch = time.date2jd(epoch_year, epoch_mon, epoch_day, epoch_hr, epoch_min, epoch_sec)[0]

        jump = (JD_i - JD_epoch) * 86400

        """Compute Constants of Perturbations"""
        raan_dot, w_dot, e_dot = PREDICT.j2dragpert(i_0, e_0, p_0, n_0, nrate_0)

        """Jump to Start Time"""
        import pdb; pdb.set_trace()
        n_s, e_s, raan_s, w_s, theta_s, M_s = PREDICT.update(jump, n_0, nrate_0, e_0, e_dot, raan_0, raan_dot, w_0, w_dot, M_0)
        i_s = i_0

        """Loop Times"""
        stp = JD_i
        while (stp < JD_f):

            """Get Satellite Info At Time"""
            r_sat_eci_turned, v_sat_eci = PREDICT.coe2rv(n_s, e_s, raan_s, w_s, theta_s, i_s, 398600.5)
            r_sat_eci = np.array([float(r_sat_eci_turned[0]), float(r_sat_eci_turned[1]), float(r_sat_eci_turned[2])])

            """Get Site Info At Time"""
            ecef2eci = COMFIX.ecef2eci(stp, lon)
            r_site_eci_turned = np.dot(ecef2eci, r_site_ecef)
            r_site_eci = np.array([r_site_eci_turned[0][0], r_site_eci_turned[1][0], r_site_eci_turned[2][0]])

            """Check Visibility At Time"""
            GST, LST = time.gstlst(stp,lon)
            import pdb; pdb.set_trace()
            visible_check = PREDICT.visible(r_sat_eci, r_site_eci, lat, lon, LST, stp)

            """Determine Where to Look At Time"""
            if (visible_check == 1):
                print("Y Script")
                rang, azm, elev = PREDICT.rhoazel(r_sat_eci, r_site_eci, lat, lon, GST)
                year, month, day, hour, minute, second = time.jd2date(stp)
                Line_out = ['Where to Look at {}/{}/{}\t {}:{}:{} JD:\n\tRange(km): {}\n\tAzmuth(deg): {}\n\tElevation(deg): {}\n'.format(int(month), int(day), int(year), int(hour), int(minute), int(second), rang, np.rad2deg(azm), np.rad2deg(elev))]
                finalfile.writelines(Line_out)

            """Update Time and Position"""
            n_s, e_s, raan_s, w_s, theta_s, M_s = PREDICT.update(time_stepsec, n_s, nrate_0, e_s, e_dot, raan_s, raan_dot, w_s, w_dot, M_s)

            stp = stp + time_stepJD
            continue
