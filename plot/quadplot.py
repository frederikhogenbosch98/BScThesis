import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt('../DGPT_DL/data.txt')


plt.plot(np.flip(x), label='phi')

plt.legend()
plt.title('phi')
plt.show()


