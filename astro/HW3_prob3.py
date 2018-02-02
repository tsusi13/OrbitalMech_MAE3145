import numpy as np
from astro import constants

def r_p(spec_mech, h, mu = 398600.5):

    ans = np.roots([2 * spec_mech, 2 * mu, -(h * h)])

    return ans
