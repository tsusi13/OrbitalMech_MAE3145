import numpy as np

mu = 1.32712e11
a1 = 1.495898e8
a2 = 4.49825e9

r1 = a1
r2 = a2
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
print("time of flight years=\n\t", tof/3600/24/365)

phi = np.pi - (np.sqrt(mu/(at*at*at))*(tof))
SPeriod = (2*np.pi)/(np.sqrt(mu/(r1*r1*r1))-np.sqrt(mu/(at*at*at)))

print("phi =\n\t", phi)
print("time hours=\n\t", SPeriod/3600/24/365)
