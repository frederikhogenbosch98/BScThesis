import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt('comparegauss.txt')
x2 = np.loadtxt('compare2cm.txt')
x4 = np.loadtxt('compare4cm.txt')
x6 = np.loadtxt('compare6cm.txt')


plt.plot(x, label='x = 0cm')
plt.plot(x2, label='x = 2cm')
plt.plot(x4, label='x = 4cm')
plt.plot(x6, label='x = 6cm', c='m')
plt.title('Flux versus energy')
plt.xlabel('Energy (MeV)')
plt.ylabel('Flux (part $cm^{-2}$$MeV^{-1}$')
plt.legend()
plt.show()
