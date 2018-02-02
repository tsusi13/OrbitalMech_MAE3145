r""" This module is designed to test the functions contained in the module COMFIX

Author: Thomas J Susi
"""

import numpy as np

from astro import COMFIX

def test_topo2rv():
    r""" Test for the TOPO2RV Function
    """

    rang = 1
    azm = 2
    elev = 3
    rang_r = 4
    azm_r = 5
    elev_r = 6
    actual_out = COMFIX.topo2rv(rang, azm, elev, rang_r, azm_r, elev_r)
    expected_out = [[-0.41198225],[-0.90019763],[0.14112001]], [[-6.501277],[-2.31079965],[-5.37547495]]

    np.testing.assert_allclose(actual_out, expected_out)

def test_lla2ecef():
    r""" Test for the LLA2ECEF Function
    """

    lat = 1
    lon = 2
    alt = 3
    actual_out = COMFIX.lla2ecef(lat, lon, alt)
    expected_out = [[-1438.17843693],[3142.47721517],[5346.2923982]]

    np.testing.assert_allclose(actual_out, expected_out)

def test_sez2ecef():
    r""" Test for the SEZ2ECEF Function
    """

    lat = 1
    lon = 2
    alt = 3
    actual_out = COMFIX.sez2ecef(lat, lon, alt)
    expected_out = [[0.84147098,0,-0.54030231],[0,1,0],[0.54030231,0,0.84147098]], [[-0.41614684,0.90929743,0],[-0.90929743,-0.41614684,0],[0,0,1]]

    np.testing.assert_allclose(actual_out, expected_out)

def test_ecef2eci():
    r""" Test for the ECEF2ECI Function
    """

    JD = 1
    lon = 2
    actual_out = COMFIX.ecef2eci(JD, lon)
    expected_out = [[-0.43320058,0.90129754,0],[-0.90129754,-0.43320058,0],[0,0,1]]

    np.testing.assert_allclose(actual_out, expected_out)
