import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')

q = np.loadtxt('test3_nhc_sho_3d_4c_q1.dat')
p = np.loadtxt('test3_nhc_sho_3d_4c_p1.dat')

plt.plot(q, p)
plt.show()