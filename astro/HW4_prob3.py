import numpy as np
from astro import RV2COE
import pdb
import matplotlib.pyplot as plt

a = 3 * 6378.14
e = 0.4
i = np.deg2rad(28.5)
raan = np.deg2rad(45)
w = np.deg2rad(90)
theta = np.deg2rad(235)
mu = 398600.5

h = np.sqrt(mu * a * (1 - (e*e)))

r = (a*(1-(e*e))) / (1+(e*np.cos(theta)))

nrg = -0.5*((mu*mu)/(h*h))*(1-e*e)

v = np.sqrt(2*(nrg+(mu/r)))

gamma = np.arctan((e*np.sin(theta))/(1+e*np.cos(theta)))

r_lvlh = np.matrix([[r], [0], [0]])

v_lvlh = np.matrix([[(mu/h)*e*np.sin(theta)], [(mu/h)*(1+e*np.cos(theta))], [0]])

E = np.arccos((a-r)/(a*e))

period = 2 * np.pi * np.sqrt((a*a*a) / mu)

change_t = np.sqrt(((a*a*a)/mu))*(E - (e*np.sin(E)))

M = E - (e*np.sin(E))

a1 = np.cos(raan)*np.cos(w) - np.sin(raan)*np.sin(w)*np.cos(i)
a2 = -np.cos(raan)*np.sin(w) - np.sin(raan)*np.cos(w)*np.cos(i)
a3 = np.sin(raan)*np.sin(i)

b4 = np.sin(raan)*np.cos(w) + np.cos(raan)*np.sin(w)*np.cos(i)
b5 = -np.sin(raan)*np.sin(w) + np.cos(raan)*np.cos(w)*np.cos(i)
b6 = -np.cos(raan)*np.sin(i)

c7 = np.sin(w)*np.sin(i)
c8 = np.cos(w)*np.sin(i)
c9 = np.cos(i)

pqw2ijk = np.matrix([[a1,a2,a3],[b4,b5,b6],[c7,c8,c9]])

lvlh2pqw = np.matrix([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])

r_pqw = np.matmul(lvlh2pqw,r_lvlh)
v_pqw = np.matmul(lvlh2pqw,v_lvlh)

r_eci = np.matmul(pqw2ijk,r_pqw)
v_eci = np.matmul(pqw2ijk,v_pqw)

print("R Vector(LVLH) =\n", r_lvlh)
print("V Vector(LVLH) =\n", v_lvlh)
print("R Vector(PQW) =\n", r_pqw)
print("V Vector(PQW) =\n", v_pqw)
print("R Vector(ECI) =\n", r_eci)
print("V Vector(ECI) =\n", v_eci)
print("R =\n\t", r)
print("V =\n\t", v)
print("Gamma =\n\t", np.rad2deg(gamma))
print("M =\n\t", 360 - np.rad2deg(M))
print("E =\n\t", 360 - np.rad2deg(E))
print("(t-T) =\n\t", (period - change_t)/3600)

RANGE = np.linspace(0, 2*np.pi, 1000)
radius = (a*(1-(e*e))) / (1+(e*np.cos(RANGE)))
plt.plot(radius*np.cos(RANGE),radius*np.sin(RANGE))
plt.plot(r*np.cos(theta),r*np.sin(theta),'bo',markersize=5)
plt.plot(0,0,'bo',markersize=20)
plt.show()
