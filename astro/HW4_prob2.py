import numpy as np
from astro import constants
from astro import RV2COE
import matplotlib.pyplot as plt

e = 1.05
r_p = 1000
mu = 398600.5

a = r_p / (e - 1)

p = abs(a * (1 - e*e))

nrg = mu / (2*a)

true_inf = np.arccos(-1 / e)

d = 2 * np.arcsin(1/e)

v_inf = abs(np.sqrt(nrg*2))

print("A\n")
print("a =\n\t", a)
print("p =\n\t", p)
print("nrg =\n\t", nrg)
print("true inf =\n\t", np.rad2deg(true_inf))
print("d =\n\t", np.rad2deg(d))
print("v inf =\n\t", v_inf)

theta_1 = np.deg2rad(90)

r_1 = (a*(e*e - 1))/(1+e*np.cos(theta_1))

v_1 = np.sqrt(2*(nrg + (mu/r_1)))

gamma = np.arctan((e * np.sin(theta_1))/(1 + (e * np.cos(theta_1))))

h = np.sqrt(mu*a*(e*e-1))

F = 2*np.arctanh(np.sqrt((e-1)/(e+1))*np.tan(theta_1/2))

M_h = e*np.sinh(F) - F

t = M_h / (((mu*mu)/(h*h*h))*((e*e - 1)**(3/2)))

r_lvlh = np.matrix([[r_1], [0], [0]])

v_lvlh = np.matrix([[(mu/h)*e*np.sin(theta_1)], [(mu/h)*(1+e*np.cos(theta_1))], [0]])

lvlh2pqw = np.matrix([[np.cos(theta_1),-np.sin(theta_1),0],[np.sin(theta_1),np.cos(theta_1),0],[0,0,1]])

r_pqw = np.matmul(lvlh2pqw,r_lvlh)

v_pqw = np.matmul(lvlh2pqw,v_lvlh)

print("B\n")
print("r1 =\n\t", r_1)
print("v1 =\n\t", v_1)
print("H1 =\n\t", np.rad2deg(F))
print("gamma =\n\t", np.rad2deg(gamma))
print("t1 =\n\t", t/60)
print("R Vector(LVLH) =\n", r_lvlh)
print("V Vector(LVLH) =\n", v_lvlh)
print("R Vector(PQW) =\n", r_pqw)
print("V Vector(PQW) =\n", v_pqw)

theta_2 = np.deg2rad(150)

F_2 = 2*np.arctanh(np.sqrt((e-1)/(e+1))*np.tan(theta_2/2))

M_h_2 = e*np.sinh(F_2) - F_2

t_2 = M_h_2 / (((mu*mu)/(h*h*h))*((e*e - 1)**(3/2)))

r_2 = (a*(e*e - 1))/(1+e*np.cos(theta_2))

print("D\n")
print("t2 =\n\t", t_2/60)
print("Time Until:\n\t", (t_2 - t)/60)
print("r2 =\n\t", r_2)

RANGE = np.linspace(np.deg2rad(-135), np.deg2rad(135))
radius = p / (1+(e*np.cos(RANGE)))
plt.axis([-5000,50000,-10000,10000])
plt.plot(radius*np.cos(RANGE),radius*np.sin(RANGE))
plt.plot(r_1+a,0,'bo',markersize=5)
plt.plot(0,0,'bo',markersize=10)
plt.plot(0,r_1,'bo',markersize=5)
plt.show()
