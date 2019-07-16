import numpy as np
import matplotlib.pyplot as plt

filename = "./test_at/test_at_middle_ho_e/test_at_middle_ho_e.chk"

times, positions, momenta = np.loadtxt(filename, unpack=True)

n, bins, patches = plt.hist(positions, bins='auto',density=True,label="calc")
plt.xlabel("Position")
plt.ylabel("Probability Density")
plt.title("Position Distribution")

x = np.linspace(min(positions), max(positions), num=1000)
y = np.exp(-4 * np.power(x, 2))/np.sqrt(0.25 * np.pi)

plt.plot(x, y, linestyle='solid', c='red', lw=3,
        alpha=0.8, label='Analytical')
plt.legend(loc='best', frameon=False)
plt.show()
