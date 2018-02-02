import numpy as np

r_earth = 6378.137
mu = 398600.5

r_B = 14000

r_A = 7000
V_A1 = 12
theta_A1 = 0
gamma_A1 = 0

a_t = .5 * (r_A + r_B)

nrg_t = -mu / (2*a_t)

v_t1 = np.sqrt(2*((mu / r_A)+(nrg_t)))

changev_a = v_t1 - V_A1

TOF = np.pi*np.sqrt(a_t**3/mu)

Period_station = 2 * np.pi * np.sqrt(r_B**3/mu)

theta_Astation = 2*np.pi*(TOF/Period_station)

time = Period_station - TOF

a_p = (mu*(time/(2*np.pi))**2)**(1/3)

r_C = r_B

nrg_p = -mu/(2*a_p)

v_t2 = np.sqrt(2*((mu / r_B)+(nrg_t)))

v_p1 = np.sqrt(2*((mu / r_C)+(nrg_p)))

changev_b1 = v_p1 - v_t2

v_C = np.sqrt(mu / r_B)

changev_b2 = v_C - v_p1

change_total = abs(changev_a) + abs(changev_b1) + abs(changev_b2)

print("Change VA: \n\t", changev_a)
print("Time A to B: \n\t", TOF)
print("Position of ISS after Hohmann: \n\t", np.rad2deg(theta_Astation))
print("Phasing Orbit Time: \n\t", time)
print("ap: \n\t", a_p)
print("RC: \n\t", r_C)
print("Change VB: \n\t", changev_b1)
print("Change VC: \n\t", changev_b2)
print("Change VT: \n\t", change_total)
