import numpy as np
from astro import RV2COE
from astro import PROPAGATE
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

e1 = 0.75
a1 = 4.5 * 6378.14
mu = 398600.5
theta1 = np.deg2rad(120)

h1 = np.sqrt(mu * a1 * (1 - (e1*e1)))

r1 = (a1*(1-(e1*e1))) / (1+(e1*np.cos(theta1)))

if (theta1 < 0):
    theta1 = 2 * np.pi + theta1
else:
    theta1 = theta1

E1 = 2 * np.arctan(np.sqrt((1-e1) / (1+e1)) * np.tan(theta1/2))

if (E1 < 0):
    E1 = 2 * np.pi + E1
else:
    E1 = E1

t_T1 = np.sqrt(((a1*a1*a1)/mu))*(E1 - (e1*np.sin(E1)))

nrg1 = -0.5*((mu*mu)/(h1*h1))*(1-e1*e1)

v1 = np.sqrt(2*(nrg1+(mu/r1)))

gamma1 = np.arctan((e1*np.sin(theta1))/(1+e1*np.cos(theta1)))

r_lvlh1 = np.matrix([[r1], [0], [0]])

v_lvlh1 = np.matrix([[(mu/h1)*e1*np.sin(theta1)], [(mu/h1)*(1+e1*np.cos(theta1))], [0]])

lvlh2pqw = np.matrix([[np.cos(theta1),-np.sin(theta1),0],[np.sin(theta1),np.cos(theta1),0],[0,0,1]])

r_pqw1 = np.matmul(lvlh2pqw,r_lvlh1)
v_pqw1 = np.matmul(lvlh2pqw,v_lvlh1)

print("r 1=\n\t", r1)
print("v 1=\n\t", v1)
print("r PQW 1=\n", r_pqw1)
print("v PQW 1=\n", v_pqw1)
print("a 1=\n\t", a1)
print("e 1=\n\t", e1)
print("gamma 1=\n\t", np.rad2deg(gamma1))
print("true anomaly 1=\n\t", np.rad2deg(theta1))
print("E 1=\n\t", np.rad2deg(E1))
print("(t-T) 1=\n\t", t_T1)

rp2 = 2.0 * 6378.14
ra2 = 6.0 * 6378.14
r2 = r1

a2 = 0.5 * (rp2 + ra2)

e2 = (ra2 - rp2) / (rp2 + ra2)

h2 = np.sqrt(mu * a2 * (1 - (e2*e2)))

theta2 = theta1

if (theta2 < 0):
    theta2 = 2 * np.pi + theta2
else:
    theta2 = theta2

E2 = 2 * np.arctan(np.sqrt((1-e2) / (1+e2)) * np.tan(theta2/2))

if (E2 < 0):
    E2 = 2 * np.pi + E2
else:
    E2 = E2

t_T2 = np.sqrt(((a2*a2*a2)/mu))*(E2 - (e2*np.sin(E2)))

nrg2 = -0.5*((mu*mu)/(h2*h2))*(1-e2*e2)

v2 = np.sqrt(2*(nrg2+(mu/r2)))

gamma2 = np.arctan((e2*np.sin(theta2))/(1+e2*np.cos(theta2)))

r_lvlh2 = np.matrix([[r2], [0], [0]])

v_lvlh2 = np.matrix([[(mu/h2)*e2*np.sin(theta2)], [(mu/h2)*(1+e2*np.cos(theta2))], [0]])

lvlh2pqw = np.matrix([[np.cos(theta2),-np.sin(theta2),0],[np.sin(theta2),np.cos(theta2),0],[0,0,1]])

r_pqw2 = np.matmul(lvlh2pqw,r_lvlh2)
v_pqw2 = np.matmul(lvlh2pqw,v_lvlh2)

print("\n\n")
print("r 2=\n\t", r2)
print("v 2=\n\t", v2)
print("r PQW 2=\n", r_pqw2)
print("v PQW 2=\n", v_pqw2)
print("a 2=\n\t", a2)
print("e 2=\n\t", e2)
print("gamma 2=\n\t", np.rad2deg(gamma2))
print("true anomaly 2=\n\t", np.rad2deg(theta2))
print("E 2=\n\t", np.rad2deg(E2))
print("(t-T) 2=\n\t", t_T2)

change_v = np.sqrt((v1*v1) + (v2*v2) - (2*v1*v2*np.cos(gamma1 - gamma2)))

beta = np.arcsin((v2 * np.sin(gamma1 - gamma2)) / change_v)

alpha = np.pi - beta

print("\n\n")
print("change v=\n\t", change_v)
print("alpha=\n\t", np.rad2deg(alpha))
print("beta=\n\t", np.rad2deg(beta))

vector = change_v * np.array([np.cos(beta)*np.cos(alpha), np.cos(beta)*np.sin(alpha), np.sin(beta)])

print("VNC \n", vector)

a3 = a2

e3 = e2

h3 = np.sqrt(mu * a3 * (1 - (e3*e3)))

theta3 = np.deg2rad(250)

if (theta3 < 0):
    theta3 = 2 * np.pi + theta3
else:
    theta3 = theta3

r3 = (a3*(1-(e3*e3))) / (1+(e3*np.cos(theta3)))

E3 = 2 * np.arctan(np.sqrt((1-e3) / (1+e3)) * np.tan(theta3/2))

if (E3 < 0):
    E3 = 2 * np.pi + E3
else:
    E3 = E3

t_T3 = np.sqrt(((a3*a3*a3)/mu))*(E3 - (e3*np.sin(E3)))

nrg3 = -0.5*((mu*mu)/(h3*h3))*(1-e3*e3)

v3 = np.sqrt(2*(nrg3+(mu/r3)))

gamma3 = np.arctan((e3*np.sin(theta3))/(1+e3*np.cos(theta3)))

r_lvlh3 = np.matrix([[r3], [0], [0]])

v_lvlh3 = np.matrix([[(mu/h3)*e3*np.sin(theta3)], [(mu/h3)*(1+e3*np.cos(theta3))], [0]])

lvlh2pqw = np.matrix([[np.cos(theta3),-np.sin(theta3),0],[np.sin(theta3),np.cos(theta3),0],[0,0,1]])

r_pqw3 = np.matmul(lvlh2pqw,r_lvlh3)
v_pqw3 = np.matmul(lvlh2pqw,v_lvlh3)

print("\n\n")
print("r 3=\n\t", r3)
print("v 3=\n\t", v3)
print("r PQW 3=\n", r_pqw3)
print("v PQW 3=\n", v_pqw3)
print("a 3=\n\t", a3)
print("e 3=\n\t", e3)
print("gamma 3=\n\t", np.rad2deg(gamma3))
print("true anomaly 3=\n\t", np.rad2deg(theta3))
print("E 3=\n\t", np.rad2deg(E3))
print("(t-T) 3=\n\t", t_T3)

time = t_T3 - t_T2

print("\n\n")
print("time taken=\n\t", time/3600)

Range = np.linspace(0, 2*np.pi, 1000)
p1 = a1*(1-(e1*e1))
p2 = a2*(1-(e2*e2))

rad1 = p1 / (1+(e1*np.cos(Range)))
rad2 = p2 / (1+(e2*np.cos(Range)))

plt.plot(rad1 * np.cos(Range), rad1 * np.sin(Range))
plt.plot(r2 * np.cos(theta2), r2 * np.sin(theta2), 'bo', markersize=5)
plt.plot(rad2 * np.cos(Range + np.deg2rad(35)), rad2 * np.sin(Range + np.deg2rad(35)))
plt.plot(0,0,'bo',markersize=20)

plt.grid()
plt.show()
