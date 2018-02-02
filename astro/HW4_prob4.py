import numpy as np
from astro import RV2COE
import matplotlib.pyplot as plt

r_earth = 6378.137

r_input = [3*r_earth, 5*r_earth, 0]
v_input = [-3.2, 2, 2.5]
mu = 398600.5

r_mag = np.linalg.norm(r_input)
v_mag = np.linalg.norm(v_input)

h_input = RV2COE.ang_momentum(r_input, v_input)

e_vector = RV2COE.eccentricity(mu, r_input, v_input, h_input)

e = np.linalg.norm(e_vector)

n_vector = RV2COE.line_of_nodes(RV2COE.unit_vector(h_input))

p = RV2COE.semi_latus_rectum(mu, h_input)

theta = RV2COE.true_anom(r_input, h_input, e_vector)

i = RV2COE.inclination(RV2COE.unit_vector(h_input))

a = RV2COE.semi_major_axis(p, np.linalg.norm(e_vector))

raan = RV2COE.R_A_A_N(n_vector)

w = RV2COE.arg_of_periapsis(n_vector, e_vector, h_input)

gamma = RV2COE.flight_ang(np.linalg.norm(e_vector), theta)

if (theta < 0):
    theta = 2 * np.pi + theta
else:
    theta = theta
if (raan < 0):
    raan = 2 * np.pi + raan
else:
    raan = raan
if (w < 0):
    w = 2 * np.pi + w
else:
    w = w
if (i < 0):
    i = 2 * np.pi + i
else:
    i = i

E = np.arccos((a-r_mag)/(a*e))

change_t = np.sqrt(((a*a*a)/mu))*(E - (e*np.sin(E)))

M = E - (e*np.sin(E))

print("r =\n\t", r_mag)
print("v =\n\t", v_mag)
print("a =\n\t", a)
print("e =\n\t", e)
print("i =\n\t", np.rad2deg(i))
print("w =\n\t", np.rad2deg(w))
print("raan =\n\t", np.rad2deg(raan))
print("gamma =\n\t", np.rad2deg(gamma))
print("true anomaly =\n\t", np.rad2deg(theta))
print("M =\n\t", np.rad2deg(M))
print("E =\n\t", np.rad2deg(E))
print("(t-T) =\n\t", change_t/60)

RANGE = np.linspace(0, 2*np.pi, 1000)
radius = p / (1+(e*np.cos(RANGE)))
plt.plot(radius*np.cos(RANGE),radius*np.sin(RANGE))
plt.plot(r_mag*np.cos(theta),r_mag*np.sin(theta),'bo',markersize=5)
plt.plot(0,0,'bo',markersize=20)
plt.show()
