import numpy as np

x = np.loadtxt('../DGPT_DL/compare.txt')
y = np.loadtxt('../DGPT_DL/data.txt')


mse = ((x[:75] - y[:75])**2).mean()
mse = ((x - y)**2).mean()

#print(x[:75])
print(mse)
