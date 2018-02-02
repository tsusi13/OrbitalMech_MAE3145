import numpy as np
import matplotlib.pyplot as plt
from astro import constants
import math

r_hat = 38208.441964
t_hat = 0
h_hat = 0
rd_hat = -0.925821
td_hat = 2.167744
hd_hat = 0

r_vector = np.array([r_hat, t_hat, h_hat])
v_vector = np.array([rd_hat, td_hat, hd_hat])

def vector_mag(vector):

    mag = np.linalg.norm(vector)

    return mag

def unit_vector(vector):

    unit_vector = vector / vector_mag(vector)

    return unit_vector

mu = 398600.5

h_vector = np.cross(r_vector, v_vector)

e_vector = np.divide(np.cross(v_vector, h_vector), mu) - np.divide(r_vector, vector_mag(r_vector))

spec_mech = ((vector_mag(v_vector) * vector_mag(v_vector)) / 2) - (mu / vector_mag(r_vector))

semi_latus = (vector_mag(h_vector) * vector_mag(h_vector)) / mu

true_ano_test = np.arccos(np.dot(unit_vector(e_vector), unit_vector(r_vector))) * 180 / np.pi
if (np.dot(r_vector, v_vector) < 0):
    true_ano = 360 - true_ano_test
else:
    true_ano = true_ano_test

semi_major = semi_latus / (1 - (vector_mag(e_vector) * vector_mag(e_vector)))

r_comp = (np.linalg.norm(e_vector) * (np.sin(true_ano))) * (mu / np.linalg.norm(h_vector))

theta_comp = (1 + (np.linalg.norm(e_vector) * (np.cos(true_ano)))) * (mu / np.linalg.norm(h_vector))

flight_ang = -np.rad2deg(np.arctan(r_comp / theta_comp))

r_p = semi_latus / (1 + vector_mag(e_vector))

r_a = semi_latus / (1 - vector_mag(e_vector))

period = 2 * np.pi * np.sqrt((semi_major * semi_major * semi_major) / mu)

"""Conversion to PQW"""
print ("PQW Vectors"),

p_hat = r_hat * np.cos(true_ano) - t_hat * np.sin(true_ano)
q_hat = r_hat * np.sin(true_ano) + t_hat * np.cos(true_ano)
w_hat = h_hat

r_pqw = np.array([p_hat, q_hat, w_hat])

print ("R: ", r_pqw),

pd_hat = rd_hat * np.cos(true_ano) - td_hat * np.sin(true_ano)
qd_hat = rd_hat * np.sin(true_ano) + td_hat * np.cos(true_ano)
wd_hat = hd_hat

v_pqw = np.array([pd_hat, qd_hat, wd_hat])

print ("V: ", v_pqw),

"""For Plotting"""
r = vector_mag(r_vector)

angle = np.linspace(0, 2 * np.pi, 1000)

rad = semi_latus / (1 + vector_mag(e_vector) * np.cos(angle))

print ("Angular Momentum: ", h_vector),
print ("e: ", e_vector),
print ("Eccentricity: ", vector_mag(e_vector)),
print ("Specific Mechanical Energy: ", spec_mech),
print ("Semi Latus Rectum: ", semi_latus),
print ("True Anomaly: ", true_ano),
print ("Semi Major Axis: ", semi_major),
print ("Flight Path Angle: ", flight_ang),
print ("Radius of Periapsis: ", r_p),
print ("Radius of Apoapsis: ", r_a),
print ("Period: ", period),
print ("Radius: ", r),

plt.plot(rad * np.cos(angle), rad * np.sin(angle))
plt.plot(r * np.cos(true_ano * np.pi/180), r * np.sin(true_ano * np.pi/180), 'bo', markersize=10)
plt.plot(-r_a, 0, 'bo', markersize=10)
plt.plot(r_p, 0, 'bo', markersize=10)
plt.plot(0, 0, 'bo', markersize=20)
plt.grid()
plt.show()
