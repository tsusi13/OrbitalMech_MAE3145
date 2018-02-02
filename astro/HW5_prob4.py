import numpy as np
from astro import RV2COE
from astro import PROPAGATE
import matplotlib.pyplot as plt

a1 = 5 * 3397.0
e1 = 0.5
i1 = np.deg2rad(30)
raan1 = np.deg2rad(45)
w1 = np.deg2rad(-60)
theta1 = np.deg2rad(120)
mu = 42828.4

changev_xyz = np.matrix([[0.1], [-0.25], [0.2]])

a1 = np.cos(raan1)*np.cos(w1) - np.sin(raan1)*np.sin(w1)*np.cos(i1)
a2 = -np.cos(raan1)*np.sin(w1) - np.sin(raan1)*np.cos(w1)*np.cos(i1)
a3 = np.sin(raan1)*np.sin(i1)

b4 = np.sin(raan1)*np.cos(w1) + np.cos(raan1)*np.sin(w1)*np.cos(i1)
b5 = -np.sin(raan1)*np.sin(w1) + np.cos(raan1)*np.cos(w1)*np.cos(i1)
b6 = -np.cos(raan1)*np.sin(i1)

c7 = np.sin(w1)*np.sin(i1)
c8 = np.cos(w1)*np.sin(i1)
c9 = np.cos(i1)

pqw2ijk = np.matrix([[a1,a2,a3],[b4,b5,b6],[c7,c8,c9]])

lvlh2pqw = np.matrix([[np.cos(theta1),-np.sin(theta1),0],[np.sin(theta1),np.cos(theta1),0],[0,0,1]])

changev_pqw  = np.matmul(pqw2ijk.transpose(),changev_xyz)
changev_lvlh = np.matmul(lvlh2pqw.transpose(),changev_pqw)

per_out = 100 * (changev_lvlh.item(2)/(abs(changev_lvlh.item(0))+abs(changev_lvlh.item(1))+abs(changev_lvlh.item(2))))

print("change v LVLH =\n", changev_lvlh)
print("Percent Out =\n", per_out)

r1v, v1v = PROPAGATE.COE2RV(a1,e1,i1,raan1,w1,theta1,mu)

print("r \n", r1v)
print("v \n", v1v)
