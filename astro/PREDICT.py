r""" This module is designed to determine the location of a satellite given radar information

Inputs Vary

Functions:
    j2dragpert
    update
    coe2rv
    visible
    rhoazel

Author: Thomas J Susi       GWU         tsusi13@gwu.edu
"""
import numpy as np
import pdb
from astro import time
from astro import planets
from kinematics import attitude

def j2dragpert(i_0, e_0, p_0, n_0, nrate_0):
    r"""Computes the rates of change that stay constant in perturbation

    Inputs:
    i_0 :(angle) Given in units of radians
    e_0 :(scalar) Given unitless
    n_0 :(angle rate) Given in units of rad/s
    nrate_0 :(angle accel) Given in units of rad/ss

    Outputs:
    raan_dot :(angle rate) Given in units of rad/s
    w_dot :(angle rate) Given in units of rad/s
    e_dot :(rate) Given in units of 1/s

    Author: Thomas J Susi       GWU         tsusi13@gwu.edu
    """

    r_earth = 6378.137

    J2 = 1.08263e-3

    n_p = n_0 * (1 + ((3/2)*J2*(r_earth / p_0)**2)*np.sqrt(1 - e_0**2)*(1 - (3/2)*np.sin(i_0)**2))

    raan_dot = (((-3/2) * J2 * (r_earth / p_0)**2) * np.cos(i_0)) * n_p
    w_dot = ((3/2) * J2 * (r_earth / p_0)**2) * (2 - (5/2)*np.sin(i_0)**2) * n_p
    e_dot = (-2/3)*((1-e_0)/n_p)*nrate_0*2

    return raan_dot, w_dot, e_dot

def update(dt, n_0, nrate_0, e_0, e_dot, raan_0, raan_dot, w_0, w_dot, M_0):
    r"""Computes the rates of change that stay constant in perturbation

    Inputs:
    dt :(scalar) Given in units of sec
    n_0 :(angle rate) Given in units of rad/s
    nrate_0 :(angle accel) Given in units of rad/ss
    e_0 :(scalar) Given unitless
    e_dot :(rate) Given in units of 1/s
    raan_0 :(angle) Given in units of radians
    raan_dot :(angle rate) Given in units of rad/s
    w_0 :(angle) Given in units of radians
    w_dot :(angle rate) Given in units of rad/s
    M_0 :(angle) Given in units of radians

    Outputs:
    n_1 :(angle rate) Given in units of rad/s
    e_1 :(scalar) Given unitless
    raan_1 :(angle) Given in units of radians
    w_1 :(angle) Given in units of radians
    theta_1 :(angle) Given in units of radians

    Author: Thomas J Susi       GWU         tsusi13@gwu.edu
    """

    n_1 = n_0 + (2*nrate_0 * dt)
    e_1 = e_0 + (e_dot * dt)
    raan_1 = raan_0 + (raan_dot * dt)
    w_1 = w_0 + (w_dot * dt)

    M_1 = M_0 + (n_0 * dt) + (nrate_0 * (dt**2))

    E_n = M_1

    check = 1

    while (check > 1e-10):

        E_1 = E_n + ((M_1 - E_n + e_1*np.sin(E_n))/(1 - e_1*np.cos(E_n)))

        check = np.absolute(E_1 - E_n)

        E_n = E_1

    theta_1 = 2 * np.arctan(np.sqrt((1 + e_1)/(1 - e_1)) * np.tan(E_1/2))

    return n_1, e_1, raan_1, w_1, theta_1, M_1

def coe2rv(n, e, raan, w, theta, i, mu):
    r"""Computes the position and velocity vectors of the satellite in ECI frame

    Inputs:
    n_1 :(angle rate) Given in units of rad/s
    e_1 :(scalar) Given unitless
    raan_1 :(angle) Given in units of radians
    w_1 :(angle) Given in units of radians
    theta_1 :(angle) Given in units of radians

    Outputs:
    r_eci :(vector) Given in ECI frame with units km
    v_eci :(vector) Given in ECI frame with units km/s

    Author: Thomas J Susi       GWU         tsusi13@gwu.edu
    """
    import pdb; pdb.set_trace()
    a = (mu / (n**2))**(1/3)

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

    return r_eci, v_eci

def visible(r_sat_eci, r_site_eci, lat, lon, GST, JD):
    r"""Checks if visibility constraints are valid for satellite information

    Inputs:
    r_sat_eci :(vector) Given in ECI frame with units km
    r_site_eci :(vector) Given in ECI frame with units km
    lat :(angle) Given in units of radians
    LST :(angle) Given in units of radians
    JD :(time) Given as Julian Date

    Outputs:
    flag :(scalar) Given as unitless, 0 = Not Visible, 1 = Visible

    Author: Thomas J Susi       GWU         tsusi13@gwu.edu
    """

    flag = 0

    r_sun = planets.sun_earth_eci(JD)[0]

    """Range and Angle"""
    p_eci = r_sat_eci - r_site_eci

    ecef2eci = ecef2eci = np.array([[np.cos(GST),-np.sin(GST),0],[np.sin(GST),np.cos(GST),0],[0,0,1]])

    p_ecef = np.dot(np.transpose(ecef2eci), p_eci)

    b2ecef = np.array([[np.cos(lon),-np.sin(lon),0],[np.sin(lon),np.cos(lon),0],[0,0,1]])

    sez2b = np.array([[np.cos(np.pi/2 - lat),0,np.sin(np.pi/2 - lat)],[0,1,0],[-np.sin(np.pi/2 - lat),0,np.cos(np.pi/2 - lat)]])

    sez2ecef = b2ecef.dot(sez2b)

    p_sez = np.dot(np.transpose(sez2ecef), p_ecef)

    rang = np.linalg.norm(p_sez)

    elev = np.arcsin(p_sez[2] / rang)

    elev_num = abs(float(elev))

    p_mag = np.linalg.norm(p_eci)

    """Darkness"""
    dark = np.dot(r_sun, r_site_eci)

    """In Sunlight"""
    sin_angle = np.dot(r_sun, r_sat_eci) / (np.linalg.norm(r_sun) * np.linalg.norm(r_sat_eci))

    cos_angle = np.dot(r_sun, r_sat_eci) / (np.linalg.norm(r_sun) * np.linalg.norm(r_sat_eci))

    angle = np.dot(r_sun, r_sat_eci) / (np.linalg.norm(r_sun) * np.linalg.norm(r_sat_eci)) + 2*np.pi

    dist = np.linalg.norm(r_sat_eci) * angle

    dist_num = float(dist)
    import pdb; pdb.set_trace()
    """If Statement"""
    if (elev_num > np.deg2rad(10) and p_mag < 1500 and dark < 0 and dist_num > 6378.137):
        flag = 1


    return flag

def rhoazel(r_sat_eci, r_site_eci, lat, lon, GST):
    r"""Computes the location of the satellite in from the site

    Inputs:
    r_sat_eci :(vector) Given in ECI frame with units km
    r_site_eci :(vector) Given in ECI frame with units km
    lat :(angle) Given in units of radians
    LST :(angle) Given in units of radians

    Outputs:
    rang :(scalar) Given in units of km
    azm :(angle) Given in units of radians
    elev :(angle) Given in units of radians

    Author: Thomas J Susi       GWU         tsusi13@gwu.edu
    """
    p_eci = r_sat_eci - r_site_eci

    ecef2eci = ecef2eci = np.array([[np.cos(GST),-np.sin(GST),0],[np.sin(GST),np.cos(GST),0],[0,0,1]])

    p_ecef = np.dot(np.transpose(ecef2eci), p_eci)

    b2ecef = np.array([[np.cos(lon),-np.sin(lon),0],[np.sin(lon),np.cos(lon),0],[0,0,1]])

    sez2b = np.array([[np.cos(np.pi/2 - lat),0,np.sin(np.pi/2 - lat)],[0,1,0],[-np.sin(np.pi/2 - lat),0,np.cos(np.pi/2 - lat)]])

    sez2ecef = b2ecef.dot(sez2b)

    p_sez = np.dot(np.transpose(sez2ecef), p_ecef)

    rang = np.linalg.norm(p_sez)

    elev = np.arcsin(p_sez[2] / rang)

    sin_azm = (p_sez[1] / (rang*np.cos(elev)))
    cos_azm = (p_sez[0] / (-rang*np.cos(elev)))
    azm = np.arctan2(sin_azm, cos_azm) + 2*np.pi

    return rang, azm, elev
