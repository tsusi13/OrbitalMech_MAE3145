r""" This module is designed to test the functions contained in the module PREDICT

Author: Thomas J Susi
"""

import numpy as np

from astro import PREDICT

def test_j2dragpert():
    r""" Test for the j2dragpert Function
    """

    i_0 = 0.5
    e_0 = 0.5
    p_0 = 10000
    n_0 = 0.001
    nrate_0 = 5e-14
    actual_out = PREDICT.j2dragpert(i_0, e_0, p_0, n_0, nrate_0)
    expected_out = -5.799757e-7, 9.420019e-7, -3.332084e-11

    np.testing.assert_allclose(actual_out, expected_out)

def test_update():
    r""" Test for the j2dragpert Function
    """

    dt = 2
    n_0 = 0.005
    nrate_0 = 0.005
    e_0 = 0.5
    e_dot = 0.1
    raan_0 = 1.2
    raan_dot = 0.2
    w_0 = 1.2
    w_dot = 0.2
    M_0 = 1
    actual_out = PREDICT.update(dt, n_0, nrate_0, e_0, e_dot, raan_0, raan_dot, w_0, w_dot, M_0)
    expected_out = .025, 0.7, 1.6, 1.6, 2.448851109, 1.03

    np.testing.assert_allclose(actual_out, expected_out)

def test_coe2rv():
    r""" Test for the j2dragpert Function
    """

    n = 0.0316227766
    e = 0.5
    raan = 1
    w = 1
    theta = 1
    i = 1
    mu = 1
    actual_out = PREDICT.coe2rv(n, e, raan, w, theta, i, mu)
    expected_out = np.matrix([[-3.7687798 ],[-0.50029911],[ 4.51804929]]), np.matrix([[-0.23816551],[-0.42423089],[-0.04485889]])

    np.testing.assert_allclose(actual_out, expected_out)

def test_visible():
    r""" Test for the j2dragpert Function
    """

    r_sat_eci = np.array([5953.919517, -171.9317157, 3358.488474])
    r_site_eci = np.array([4931.448597, -606.6268224, 3885.876923])
    lat = 0.679369
    lon = -1.344898
    GST = 2.0713023107419848
    JD = 2458088.634722
    actual_out = PREDICT.visible(r_sat_eci, r_site_eci, lat, lon, GST, JD)
    expected_out = 1

    np.testing.assert_allclose(actual_out, expected_out)
