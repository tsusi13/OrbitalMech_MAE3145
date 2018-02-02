r""" This module is designed to determine a multitude of data points for a given orbit

Inputs Vary

Functions:
    Update
    COE2RV

Author: Thomas J Susi       GWU         tsusi13@gwu.edu
"""
import numpy as np
import pdb

def update(a_i, e_i, i_i, raan_i, w_i, theta_i, change_t, mu):
    r"""Computes the updated COEs given the time change

    Inputs:
    a_i :(scalar) Given in units of km
    e_i :(scalar) Given unitless
    i_i :(angle) Given in units of radians
    raan_i :(angle) Given in units of radians
    w_i :(angle) Given in units of radians
    theta_i :(angle) Given in units of radians
    change_t :(scalar) Givens in units of seconds
    mu :(scalar) Given in units kmkmkm/ss

    Outputs:
    a_f :(scalar) Given in units of km
    e_f :(scalar) Given unitless
    i_f :(angle) Given in units of radians
    raan_f :(angle) Given in units of radians
    w_f :(angle) Given in units of radians
    theta_f :(angle) Given in units of radians

    Author: Thomas J Susi       GWU         tsusi13@gwu.edu
    """

    E_i = 2 * np.arctan(np.sqrt((1-e_i) / (1+e_i)) * np.tan(theta_i/2))

    if (E_i < 0):
        E_i = 2 * np.pi + E_i
    else:
        E_i = E_i

    t_T_i = np.sqrt(((a_i*a_i*a_i)/mu))*(E_i - (e_i*np.sin(E_i)))

    M_i = E_i - (e_i*np.sin(E_i))

    if (M_i < 0):
        M_i = 2 * np.pi + M_i
    else:
        M_i = M_i

    print("Initial Time Since Periapsis in Seconds:\n\t", str(t_T_i), "\n")
    print("Initial Eccentric Anomaly in Degrees:\n\t", str(np.rad2deg(E_i)), "\n")
    print("Initial Mean Anomaly in Degrees:\n\t", str(np.rad2deg(M_i)), "\n")

    M_f = (np.sqrt(mu/(a_i*a_i*a_i))*change_t) + M_i

    r"""Newtons Method to determine final Eccentric Anomaly"""

    E_n = M_f

    check = 1

    while (check > 1e-10):

        E_f = E_n + ((M_f - E_n + e_i*np.sin(E_n))/(1 - e_i*np.cos(E_n)))

        check = np.absolute(E_f - E_n)

        E_n = E_f

    t_T_f = np.sqrt(((a_i*a_i*a_i)/mu))*(E_f - (e_i*np.sin(E_f)))

    theta_f = 2 * np.arctan(np.sqrt((1 + e_i)/(1 - e_i)) * np.tan(E_f/2))

    if (theta_f < 0):
        theta_f = 2 * np.pi + theta_f
    else:
        theta_f = theta_f

    a_f = a_i

    e_f = e_i

    i_f = i_i

    w_f = w_i

    raan_f = raan_i

    print("Final Time Since Periapsis in Seconds:\n\t", str(t_T_f), "\n")
    print("Final Eccentric Anomaly in Degrees:\n\t", str(np.rad2deg(E_f)), "\n")
    print("Final Mean Anomaly in Degrees:\n\t", str(np.rad2deg(M_f)), "\n")
    print("Final Semi Major Axis in Kilometers:\n\t", str(a_f), "\n")
    print("Final Eccentricity:\n\t", str(e_f), "\n")
    print("Final Inclination in Degrees:\n\t", str(np.rad2deg(i_f)), "\n")
    print("Final RAAN in Degrees:\n\t", str(np.rad2deg(raan_f)), "\n")
    print("Final Argument of Periapsis in Degrees:\n\t", str(np.rad2deg(w_f)), "\n")
    print("Final True Anomaly in Degrees:\n\t", str(np.rad2deg(theta_f)), "\n")

    return a_f, e_f, i_f, raan_f, w_f, theta_f

def COE2RV(a, e, i, raan, w, theta, mu):
    r"""Computes the radius and velocity vectors given COEs

    Inputs:
    a :(scalar) Given in units of km
    e :(scalar) Given unitless
    i :(angle) Given in units of radians
    raan :(angle) Given in units of radians
    w :(angle) Given in units of radians
    theta :(angle) Given in units of radians
    mu :(scalar) Given in units kmkmkm/ss

    Outputs:
    r_lvlh :(vector) in LVLH with units in km
    v_lvlh :(vector) in LVLH with units in km/s
    r_pqw :(vector) in PQW with units in km
    v_pqw :(vector) in PQW with units in km/s
    r_eci :(vector) in ECI with units in km
    v_eci :(vector) in ECI with units in km/s

    Author: Thomas J Susi       GWU         tsusi13@gwu.edu
    """

    h = np.sqrt(mu * a * (1 - (e*e)))

    r = (a*(1-(e*e))) / (1+(e*np.cos(theta)))

    nrg = -0.5*((mu*mu)/(h*h))*(1-e*e)

    v = np.sqrt(2*(nrg+(mu/r)))

    gamma = np.arctan((e*np.sin(theta))/(1+e*np.cos(theta)))

    r_lvlh = np.matrix([[r], [0], [0]])

    v_lvlh = np.matrix([[(mu/h)*e*np.sin(theta)], [(mu/h)*(1+e*np.cos(theta))], [0]])

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

    return r_eci, v_eci
