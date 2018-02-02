import numpy as np
import matplotlib.pyplot as plt
from astro import constants
"""
print ("Give the IJK values for the R vector: "),
r_vector = np.array([input(), input(),input()])

print ("Give the IJK values for the V vector: "),
v_vector = np.array([input(), input(),input()])
"""

r_vector = np.array([-7546026.6144396, -3717105.21901527, -1515557.34280287])
v_vector = np.array([0.89506649, -0.33312074, 0.21519571])

def vector_mag(vector):

    mag = np.linalg.norm(vector)

    return mag

def unit_vector(vector):

    unit_vector = vector / vector_mag(vector)

    return unit_vector

mu = 3.79313e7

h_vector = np.cross(r_vector, v_vector)

e_vector = ((1 / mu) * ((np.linalg.norm(v_vector)**2 - (mu / np.linalg.norm(r_vector))) * r_vector - (r_vector.dot(v_vector) * v_vector)))

spec_mech = ((vector_mag(v_vector) * vector_mag(v_vector)) / 2) - (mu / vector_mag(r_vector))

semi_latus = (vector_mag(h_vector) * vector_mag(h_vector)) / mu

line_nodes = np.cross([0,0,1], unit_vector(h_vector))

inclin = np.arccos(np.dot([0,0,1], unit_vector(h_vector))) * 180 / np.pi

true_ano = np.arccos(np.dot(unit_vector(e_vector), unit_vector(r_vector))) * 180 / np.pi

semi_major = semi_latus / (1 - (vector_mag(e_vector) * vector_mag(e_vector)))

RAAN = np.arccos(np.dot([1,0,0], line_nodes)) * 180 / np.pi

arg_peri = np.arccos(np.dot(line_nodes, unit_vector(e_vector))) * 180 / np.pi

flight_ang = np.arccos(h_vector / np.dot(r_vector, v_vector))

r_p = semi_latus / (1 + vector_mag(e_vector))

r_a = semi_latus / (1 - vector_mag(e_vector))

period = 2 * np.pi * np.sqrt((semi_major * semi_major * semi_major) / mu)

"""For Plotting"""
r = vector_mag(r_vector)

angle = np.linspace(0, 2 * np.pi, 1000)

rad = semi_latus / (1 + vector_mag(e_vector) * np.cos(angle))

print ("Angular Momentum: ", h_vector),
print ("e: ", e_vector),
print ("Eccentricity: ", vector_mag(e_vector)),
print ("Specific Mechanical Energy: ", spec_mech),
print ("Semi Latus Rectum: ", semi_latus),
print ("Line of Nodes: ", line_nodes),
print ("Inclination: ", inclin),
print ("True Anomaly: ", true_ano),
print ("Semi Major Axis: ", semi_major),
print ("RAAN: ", RAAN),
print ("Argument of Periapsis: ", arg_peri),
print ("Flight Path Angle: ", flight_ang),
print ("Radius of Periapsis: ", r_p),
print ("Radius of Apoapsis: ", r_a),
print ("Period: ", period),
print ("Radius: ", r),

"""Plotting the Orbits in both Systems"""

fig, ax = plt.subplots()
ax.plot(rad * np.cos(angle), rad * np.sin(angle))
ax.plot(r * np.cos(true_ano * np.pi/180), r * np.sin(true_ano * np.pi/180), 'ro', markersize = 10)
ax.plot(0, 0, 'bo', markersize = 20)
ax.plot(r_p, 0 , 'go', markersize = 5)
ax.set_xlabel('$\hat p$')
ax.set_ylabel('$\hat q$')
plt.grid()
plt.show()
