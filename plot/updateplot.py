import numpy as np
import matplotlib.pyplot as plt

before = np.loadtxt('../DGPT_DL/beforeupdate.txt')
after = np.loadtxt('../DGPT_DL/afterupdate.txt')

plt.plot(before, label='before update')
plt.plot(after, label='after update')
plt.legend()
plt.title('Subset updating')
plt.xlabel('Energy groups')
plt.ylabel('Flux (part $cm^{-2}MeV^{-1}$')
plt.show()
