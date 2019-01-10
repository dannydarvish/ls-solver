import matplotlib.pyplot as plt
import numpy as np
import re

energy_list = []
delta_list = []
# for E in range(300,900,10) + range(900,1000,1) + range(1000,1500,10):
# for E in range(1000,1201,5):
# for E in range(750,1201,10):
for E in range(300,1510,10):
    if E == 700:
        E = 690
    try:
        f = open('jobs/E_%d.log' % E)
    except IOError:
        print('E = %f file not present' % E)
        continue
    energy_list.append(E)
    lines = f.read().splitlines()
    for line in lines:
        if 'delta' in line:
            delta = float(re.search('delta = (\d*.\d*).', line).group(1))
            # if E < 1010:
            #     delta += np.pi
            # if E == 1000:
            #     delta += np.pi
            if E >= 1000:
                delta += np.pi
            delta_list.append(delta*180.0/np.pi)

plt.plot(energy_list, delta_list)
plt.xlabel(r'$E$ (MeV)')
plt.ylabel(r'Re$\delta_{\pi\pi}$ (degrees)')
plt.show()