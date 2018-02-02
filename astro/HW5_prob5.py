import numpy as np
from astro import RV2COE
from astro import PROPAGATE
import matplotlib.pyplot as plt

mu = 398600.5
r1 = 1.25 * 6378.14
r2 = 6.6 * 6378.14

a1 = r1
e1 = 0
gamma1 = 0

v1 = np.sqrt(2 * ((mu / r1) - (mu / (2 * a1))))

at = .5 * (r1 + r2)
et = 1 - (r1 / at)

vt1 = np.sqrt(2 * ((mu / r1) - (mu / (2 * at))))
gammat1 = 0

changev1 = vt1 - v1

vt2 = np.sqrt(2 * ((mu / r2) - (mu / (2 * at))))
gammat2 = 0

v2 = np.sqrt(mu / r2)

changev2 = v2 - vt2

changev = changev2 + changev1

tof = np.pi * np.sqrt((at*at*at)/mu)

print("change 1=\n\t", changev1)
print("change 2=\n\t", changev2)
print("change total=\n\t", changev)
print("time of flight=\n\t", tof/3600)

phi = np.pi - (np.sqrt(mu/(at*at*at))*(tof))
time = (2*np.pi)/(np.sqrt(mu/(r1*r1*r1))-np.sqrt(mu/(at*at*at)))

print("phi =\n\t", phi)
print("time =\n\t", time/3600)
