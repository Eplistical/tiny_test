#!/usr/bin/env python3

import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns


fig = plt.figure()
ax = fig.gca()

a = np.loadtxt('1')
#plt.plot(a[:,1], a[:,4], label='x')
plt.plot(a[:,1], a[:,2], label='xy-plane')

plt.legend()
plt.show()
