import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



x = np.loadtxt('../DGPT_DL/data.txt')
y = np.loadtxt('../DGPT_DL/dataphi.txt')

print(len(x))
print(len(y))

print(y)

plt.plot(x, label='before update')
plt.plot(y, label='after update')
plt.legend()
plt.title('phi')
plt.show()

