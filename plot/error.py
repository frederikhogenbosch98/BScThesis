import numpy as np

x = np.loadtxt('comparecoef.txt')
y = np.loadtxt('../DGPT_DL/boundedcoef.txt')



mse = ((x - y)**2).mean()
pe2 = mse/float((x**2).mean())
sqe = np.sqrt(pe2)
#print(x[:75])
print('mean squared error (mse): ',mse)
print('norm mse (nmse): ', pe2)
print('sq. r. of nmse (sqnmse): ', sqe)
