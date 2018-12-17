import matplotlib.pyplot as plt
import numpy as np

results = np.loadtxt('delta_pipi_re_1b1c.txt')
plt.plot(results[:,0],results[:,1])
plt.ylim([0,400])
plt.xlim([300,1000])
plt.show()
