set xrange[50:150]
#set yrange[0.81:0.825]
T = 0.1/2
sigma0=10
t0 = sigma0**2/(4*T)
t = t0+1000
sigma = sqrt(4*T*t)
plot 'magweg' w l,(sigma0/sigma)*exp(-(x-100)**2/sigma**2)
pause -1
