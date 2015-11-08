# test suite for Pade_multi.py

from Pade_multi import Pade_fit, Pade_eval
import numpy as np
# let's see how well it does with a polynomial
X1 = np.linspace(-1.0,1.0, 20)
X2 = np.copy(X1)
X1v,X2v = np.meshgrid(X1,X2)
Yv = X1v * X1v + X2v * X2v
X1 = np.ravel(X1v)
X2 = np.ravel(X2v)
Y = np.ravel(Yv)

XV = np.vstack((X1,X2,Y))
XV = XV.T
ans, ae, be = Pade_fit(XV,2,2)


YA = []
for i in range(XV.shape[0])  :
    YA.append(Pade_eval( XV[i,:-1], ans, ae, be ))
YA = np.array(YA)
#YA = YA.reshape(X1v.shape)


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X1,X2,Y,c=u'r')
ax.scatter(X1,X2,YA)
plt.show()
