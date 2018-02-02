import numpy as np

mu = 4902.8
e = 0
a1 = r1 = 1737.5 + 100

h1 = np.sqrt(mu * a1 * (1 - (e*e)))

nrg1 = -0.5*((mu*mu)/(h1*h1))*(1-e*e)

v1 = np.sqrt(2*(nrg1+(mu/r1)))

change_v1 = 2 * v1 * np.sin((np.deg2rad(90))/2)

print("Change in V1:\n\t", change_v1)

"""Hohmann 1"""
r2 = a2 = 1737.5 + 17000

ata = .5 * (r1 + r2)

nrga = -mu / (2*ata)
vt1a = np.sqrt(2 * (mu / r1 + nrga))
gammat1a = 0

changev1a = vt1a - v1

vt2a = np.sqrt(2 * (mu / r2 + nrga))
gammat2a = 0

v2a = np.sqrt(mu / r2)

changev2a = v2a - vt2a

change_v2a = changev1a

tofa = np.pi * np.sqrt((ata*ata*ata)/mu)

"""Plane Change"""

change_v2b = 2 * vt2a * np.sin((np.deg2rad(45)))

"""Hohmann 2"""
atc = .5 * (r2 + r1)
etc = 1 - (r1 / atc)

nrgc = nrga

v2 = np.sqrt((2*(mu / r2 + nrgc)))

vt1c = np.sqrt(2 * (mu / r1 +nrgc))
gammat1c = 0

changev1c = v2 - vt1c

change_v2c = vt1c - v1

tofc = np.pi * np.sqrt((ata*ata*ata)/mu)

"""Addition"""
change_v2 = change_v2a + change_v2b + change_v2c

toft = tofa + tofc

print("Change of V2:\n\t", change_v2)
print("TOF Bi-Elliptic Days:\n\t", toft/3600/24)
