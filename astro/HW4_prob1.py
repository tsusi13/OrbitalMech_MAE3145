import numpy as np
from astro import constants
from astro import RV2COE

e = 0.9671429
a_ua = 17.834144
true_anom = np.deg2rad(float(input("True Anomaly?: ")))

mu = 1.32712e11

ua = 1.495898e8

a = a_ua * ua

h = np.sqrt(mu * a * (1 - (e*e)))

r_p = a * (1 - e)

r_a = a * (1 + e)

period = 2 * np.pi * np.sqrt((a*a*a) / mu)

nrg = - (mu / (2 * a))

r = (a*(1-(e*e))) / (1+(e*np.cos(true_anom)))

v = np.sqrt(2*(nrg + (mu/r)))

flight_ang = np.arctan((e * np.sin(true_anom))/(1 + (e * np.cos(true_anom))))

E = np.arccos((a-r)/(a*e))

change_t = np.sqrt(((a*a*a)/mu))*(E - (e*np.sin(E)))

print("r =\n\t", r/ua)
print("v =\n\t", v)
print("gamma =\n\t", np.rad2deg(flight_ang))
print("energy =\n\t", nrg)
print("h =\n\t", h)
print("period =\n\t", period/31536000)
print("r_a =\n\t", r_a/ua)
print("r_p =\n\t", r_p/ua)
print("E =\n\t", 360 - np.rad2deg(E))
print("(t-T) =\n\t", period/31536000 - change_t/31536000)
