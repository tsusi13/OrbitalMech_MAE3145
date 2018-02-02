import numpy as np
from astro import RV2COE
from astro import PROPAGATE
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

mu = 398600.5
r_1 = 4 * 6378.14
v_1 = 4.54
gamma_1 = np.deg2rad(-40)

nrg_1 = ((v_1 * v_1) / 2) - (mu / r_1)

a_1 = -mu / (2 * nrg_1)

e_1 = np.sqrt(-(((r_1 * v_1 * np.cos(gamma_1))**2) / (mu * a_1)) + 1)

period_1 = 2 * np.pi * np.sqrt((a_1**3) / mu)

r_p_1 = a_1 * (1 - e_1)

r_a_1 = a_1 * (1 + e_1)

theta_1 = -np.arccos(((((1 - (e_1**2)) * a_1) / r_1) - 1) / e_1)

if (theta_1 < 0):
    theta_1 = 2 * np.pi + theta_1
else:
    theta_1 = theta_1

E_1 = 2 * np.arctan(np.sqrt((1-e_1) / (1+e_1)) * np.tan(theta_1/2))

if (E_1 < 0):
    E_1 = 2 * np.pi + E_1
else:
    E_1 = E_1

t_T_1 = np.sqrt(((a_1*a_1*a_1)/mu))*(E_1 - (e_1*np.sin(E_1)))

p_1 = a_1*(1-(e_1*e_1))

print("r =\n\t", r_1)
print("v =\n\t", v_1)
print("a =\n\t", a_1)
print("e =\n\t", e_1)
print("rp =\n\t", r_p_1)
print("ra =\n\t", r_a_1)
print("NRG =\n\t", nrg_1)
print("period =\n\t", period_1)
print("gamma =\n\t", np.rad2deg(gamma_1))
print("true anomaly =\n\t", np.rad2deg(theta_1))
print("E =\n\t", np.rad2deg(E_1))
print("(t-T) =\n\t", t_T_1)

a_2, e_2, i_2, raan_2, w_2, theta_2 = PROPAGATE.update(a_1, e_1, 1, 1, 1, theta_1, (8.5 * 3600), mu)

r_2 = (a_2 * (1 - e_2**2)) / (1 + e_2*np.cos(theta_2))

v_2 = np.sqrt(2 * ((mu / r_2) - (mu / (2 * a_2))))

gamma_2 = np.arccos(np.sqrt(mu * a_2 * (1 - e_2**2)) / (r_2 * v_2))

E_2 = 2 * np.arctan(np.sqrt((1-e_2) / (1+e_2)) * np.tan(theta_2/2))

if (E_2 < 0):
    E_2 = 2 * np.pi + E_2
else:
    E_2 = E_2

t_T_2 = np.sqrt(((a_2*a_2*a_2)/mu))*(E_2 - (e_2*np.sin(E_2)))

print("a 2 =\n\t", a_2)
print("e 2 =\n\t", e_2)
print("true anomaly 2 =\n\t", np.rad2deg(theta_2))
print("E 2 =\n\t", np.rad2deg(E_2))
print("(t-T) 2 =\n\t", t_T_2)
print("r 2 =\n\t", r_2)
print("v 2 =\n\t", v_2)
print("Gamma 2 =\n\t", np.rad2deg(gamma_2))

alpha = np.deg2rad(30)

change_v = 1.2

r_f = r_2

v_f = np.sqrt((v_2**2) + (change_v**2) - (2 * v_2 * change_v * np.cos(np.pi - alpha)))

gamma_f = np.arcsin((change_v / v_f) * np.sin(np.deg2rad(-150)))

nrg_f = ((v_f * v_f) / 2) - (mu / r_f)

a_f = -mu / (2 * nrg_f)

e_f = np.sqrt(-(((r_2 * v_f * np.cos(gamma_f))**2) / (mu * a_f)) + 1)

period_f = 2 * np.pi * np.sqrt((a_f**3) / mu)

r_p_f = a_f * (1 - e_f)

r_a_f = a_f * (1 + e_f)

theta_f = -np.arccos(((((1 - (e_f**2)) * a_f) / r_f) - 1) / e_f)

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

print("r final=\n\t", r_2)
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

deltw = theta_2 - theta_f
print("change w =\n\t", np.rad2deg(deltw))

Range = np.linspace(0, 2*np.pi, 1000)
rad1 = p_1 / (1+(e_1*np.cos(Range)))
rad2 = p_f / (1+(e_f*np.cos(Range)))

plt.plot(rad1 * np.cos(Range), rad1 * np.sin(Range))
plt.plot(r_2 * np.cos(theta_2), r_2 * np.sin(theta_2), 'bo', markersize=5)
plt.plot(rad2 * np.cos(Range + deltw), rad2 * np.sin(Range + deltw))
plt.plot(0,0,'bo',markersize=20)
plt.plot(r_p_1, 0, 'bo', markersize=3)
plt.plot(-r_a_1, 0, 'bo', markersize=3)
plt.grid()
plt.show()
