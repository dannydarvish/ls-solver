import matplotlib.pyplot as plt
import numpy as np
import re

energy_list = []
delta_list = []
# for E in range(300,1600,100):
# for E in [300,400,500,600,690,800,900,1000,1200,1300,1400,1500]:
for E in range(300,1010,10):
    if E == 700:
        E = 690
    energy_list.append(E)
    f = open('jobs/E_%d.log' % E)
    lines = f.read().splitlines()
    for line in lines:
        if 'delta' in line:
            delta = float(re.search('delta = (\d*.\d*).', line).group(1))
    delta_list.append(delta*180.0/np.pi)

plt.plot(energy_list, delta_list)
plt.show()