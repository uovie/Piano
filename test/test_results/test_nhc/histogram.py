import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')


def func(x):
    return  2 * x * x * np.exp(-x * x / 2)/np.sqrt(2 * np.pi)


d = np.loadtxt('E:/Users/Leoy/Desktop/test_pimd_sho_3d_q1.dat')

#hist, bin_edges = np.histogram(d)

#  matplotlib.axes.Axes.hist() 方法的接口

fig, ax = plt.subplots()
n, bins, patches = ax.hist(x=d, bins='auto',density=True,label="calc")
#color='#0504aa',alpha=0.7, rwidth=0.85)
#plt.grid(axis='y', alpha=0.75)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('My Very Own Histogram')
plt.text(6, 0.1, r'$\mu=0, \sigma=1$')
maxfreq = n.max()
x = sorted(d)
y = []
for xi in x:
    y.append(func(xi))

ax.plot(x, y, linestyle='solid', c='red', lw=3,
        alpha=0.8, label='Analytical')
ax.legend(loc='best', frameon=False)
plt.show()


