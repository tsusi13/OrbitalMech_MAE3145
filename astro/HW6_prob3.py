import numpy as np
import matplotlib.pyplot as plt

e = 0.4
a = 6 * 6378.14
mu = 398600.5

theta = np.deg2rad(90)

alpha = np.deg2rad(-60)
beta = np.pi - alpha
changevmag = 0.75

gamma = np.arctan((e*np.sin(theta))/(1+e*np.cos(theta)))

changev_vnc = changevmag * np.array([np.cos(beta)*np.cos(alpha), np.cos(beta)*np.sin(alpha), np.sin(beta)])

vnc2lvlh = np.array([[np.sin(gamma), np.cos(gamma), 0],[np.cos(gamma), -np.sin(gamma), 0],[0,0,1]])

changev_lvlh = np.dot(vnc2lvlh,changev_vnc)

lvlh2pqw = np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])

changev_pqw = np.dot(lvlh2pqw,changev_lvlh)

print("Change V VNC: \n\t", changev_vnc)
print("Change V LVLH: \n\t", changev_lvlh)
print("Change V PQW: \n\t", changev_pqw)

r_f = (a*(1-e**2))/(1 + e*np.cos(theta))
v1 = np.sqrt(2*((mu/r_f)-(mu/(2*a))))

h1 = r_f*v1*np.cos(gamma)

v_f = np.sqrt((v1**2) + (changevmag**2) - (2*v1*changevmag*np.cos(np.pi - alpha)))

nrg_f = ((v_f * v_f) / 2) - (mu / r_f)

a_f = -mu / (2 * nrg_f)

gamma_f = np.arccos((changevmag**2 - v_f**2 - v1**2) / (-2 * v_f * v1)) - gamma

h_f = r_f*v_f*np.cos(gamma_f)

p_f = (h_f**2)/mu

e_f = np.sqrt(-(((r_f * v_f * np.cos(gamma_f))**2) / (mu * a_f)) + 1)

period_f = 2 * np.pi * np.sqrt((a_f**3) / mu)

r_p_f = a_f * (1 - e_f)

r_a_f = a_f * (1 + e_f)

theta_f = np.arccos(((((1 - (e_f**2)) * a_f) / r_f) - 1) / e_f)

if (theta_f < 0):
    theta_f = 2 * np.pi + theta_f
else:
    theta_f = theta_f

E_f = 2 * np.arctan(np.sqrt((1-e_f) / (1+e_f)) * np.tan(theta_f/2))

if (E_f < 0):
    E_f = 2 * np.pi + E_f
else:
    E_f = E_f

t_T_f = np.sqrt(((a_f*a_f*a_f)/mu))*(E_f - (e_f*np.sin(E_f)))

p_f = a_f*(1-(e_f*e_f))

change_w = theta - theta_f

print("Change @ 90: \n")
print("r final=\n\t", r_f)
print("v final=\n\t", v_f)
print("a final=\n\t", a_f)
print("e final=\n\t", e_f)
print("rp final=\n\t", r_p_f)
print("ra final=\n\t", r_a_f)
print("NRG final=\n\t", nrg_f)
print("period final=\n\t", period_f)
print("gamma final=\n\t", np.rad2deg(gamma_f))
print("true anomaly final=\n\t", np.rad2deg(theta_f))
print("E final=\n\t", np.rad2deg(E_f))
print("(t-T) final=\n\t", t_T_f)
print("Change W=\n\t", np.rad2deg(change_w))

"""Wait Time to 270"""
theta_90 = np.deg2rad(90)
E_90 = 2 * np.arctan(np.sqrt((1-e) / (1+e)) * np.tan(theta_90/2))

if (E_90 < 0):
    E_90 = 2 * np.pi + E_90
else:
    E_90 = E_90

t_T_90 = np.sqrt(((a*a*a)/mu))*(E_90 - (e*np.sin(E_90)))

E_270 = np.deg2rad(270)

theta_270 = np.arccos((np.cos(E_270) - e)/(1 - e*np.cos(E_270)))

if (E_270 < 0):
    E_270 = 2 * np.pi + E_270
else:
    E_270 = E_270

t_T_270 = np.sqrt(((a*a*a)/mu))*(E_270 - (e*np.sin(E_270)))

time = t_T_270 - t_T_90

print("Time from 90 to 270: \n\t", time)

r_f_2 = (a*(1-e**2))/(1 + e*np.cos(theta_270))
v1_2 = np.sqrt(2*((mu/r_f)-(mu/(2*a))))

v_f_2 = np.sqrt((v1_2**2) + (changevmag**2) - (2*v1_2*changevmag*np.cos(np.pi - alpha)))

nrg_f_2 = ((v_f_2 * v_f_2) / 2) - (mu / r_f_2)

a_f_2 = -mu / (2 * nrg_f_2)

gamma_2 = np.arctan((e*np.sin(theta_270))/(1+e*np.cos(theta_270)))

gamma_f_2 = np.arccos((changevmag**2 - v_f_2**2 - v1_2**2) / (-2 * v_f_2 * v1_2)) - gamma_2

h_f_2 = r_f_2*v_f_2*np.cos(gamma_f_2)

p_f_2 = (h_f_2**2)/mu

e_f_2 = np.sqrt(-(((r_f_2 * v_f_2 * np.cos(gamma_f_2))**2) / (mu * a_f_2)) + 1)

period_f_2 = 2 * np.pi * np.sqrt((a_f_2**3) / mu)

r_p_f_2 = a_f_2 * (1 - e_f_2)

r_a_f_2 = a_f_2 * (1 + e_f_2)

theta_f_2 = np.arccos(((((1 - (e_f_2**2)) * a_f_2) / r_f_2) - 1) / e_f_2)

if (theta_f_2 < 0):
    theta_f_2 = 2 * np.pi + theta_f_2
else:
    theta_f_2 = theta_f_2

E_f_2 = 2 * np.arctan(np.sqrt((1-e_f_2) / (1+e_f_2)) * np.tan(theta_f_2/2))

if (E_f_2 < 0):
    E_f_2 = 2 * np.pi + E_f_2
else:
    E_f_2 = E_f_2

t_T_f_2 = np.sqrt(((a_f_2*a_f_2*a_f_2)/mu))*(E_f_2 - (e_f_2*np.sin(E_f_2)))

p_f_2 = a_f_2*(1-(e_f_2*e_f_2))

change_w_2 = theta_270 - theta_f_2

print("Change @ 270: \n")
print("r final=\n\t", r_f_2)
print("v final=\n\t", v_f_2)
print("a final=\n\t", a_f_2)
print("e final=\n\t", e_f_2)
print("rp final=\n\t", r_p_f_2)
print("ra final=\n\t", r_a_f_2)
print("NRG final=\n\t", nrg_f_2)
print("period final=\n\t", period_f_2)
print("gamma final=\n\t", np.rad2deg(gamma_f_2))
print("true anomaly final=\n\t", np.rad2deg(theta_f_2))
print("E final=\n\t", np.rad2deg(E_f_2))
print("(t-T) final=\n\t", t_T_f_2)
print("Change W=\n\t", np.rad2deg(change_w_2))

"""Plot"""
p_1 = (h1*h1)/mu
Range = np.linspace(0,2*np.pi,1000)
rad1 = p_1/(1+(e*np.cos(Range)))
plt.plot(rad1*np.cos(Range), rad1*np.sin(Range))
plt.plot(r_f*np.cos(theta), r_f*np.sin(theta), 'bo', markersize=5)
plt.plot(0,0,'bo', markersize=15)

radf = p_f/(1+(e_f*np.cos(Range)))
plt.plot(radf*np.cos(Range + change_w), radf*np.sin(Range + change_w))

radf_2 = p_f_2/(1+(e_f_2*np.cos(Range)))
plt.plot(radf_2*np.cos(Range+change_w_2), radf_2*np.sin(Range+change_w_2))

plt.grid()
plt.show()
