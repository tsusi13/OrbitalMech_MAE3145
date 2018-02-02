r""" This module is designed to test the functions contained in the module PROPAGATE

Author: Thomas J Susi
"""

import numpy as np

from astro import PROPAGATE

def test_update():
    r""" Test for the Update Function
    """

    a = 15300
    e = 0.37254902
    i = 1
    raan = 1
    w = 1
    theta = np.deg2rad(120)
    change_t = 5340
    mu = 398600.5
    actual_out = PROPAGATE.update(a, e, i, raan, w, theta, change_t, mu)
    expected_out = 15300, 0.37254902, 1, 1, 1, 3.141580

    np.testing.assert_allclose(actual_out, expected_out)

def test_COE2RV():
    r""" Test for the COE2RV Function
    """

    a = -0.0672699547601
    e = 12.7101202659
    i = np.deg2rad(114.094842552)
    raan = np.deg2rad(206.565051177)
    w = np.deg2rad(37.0961328341)
    theta = np.deg2rad(81.4646923837)
    mu = 5
    actual_out = PROPAGATE.COE2RV(a, e, i, raan, w, theta, mu)
    expected_out = [1,2,3], [4,5,6]
