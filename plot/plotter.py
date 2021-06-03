import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



x = np.loadtxt('../DGPT_DL/data.txt')
y = np.loadtxt('../DGPT_DL/dataphi.txt')

print(len(x))
print(len(y))

print(y)

#plt.plot(np.flip(x), label='end')
for k in range(len(y)):
    if k%1==0:
        plt.plot(y[k], label='y'+str(k))


#plt.plot(y, label='y')
#plt.plot(y[1], label='y1')
plt.legend()
plt.title('phi')
plt.show()

