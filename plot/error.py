import numpy as np

x = np.loadtxt('compareclean.txt')
y = np.loadtxt('../DGPT_DL/data.txt')


mse = ((x[:75] - y[:75])**2).mean()
mse = ((x - y)**2).mean()
pe = mse/float(x.size**2)
pe2 = mse/float((x**2).mean())
#print(x[:75])
print('mean squared error: ',mse)
print('percentage error: ',pe)
print('diff error: ', pe2)
