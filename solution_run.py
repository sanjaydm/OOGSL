x = [1.03218082032235,1.24604460897023,1.33909468691213,1.35263245425559,1.3177913304116,1.25716426513867,1.18633831921892,1.11542260840888,1.05043416468788,0.994461434292521,0.948593583222092,0.912637792554866,0.885656614580798,0.866357319600424,0.853362331534024,0.845386561030597,0.841343977851758,0.84040200444613,0.841998279738827,0.845830166083338,0.851823253739638,0.860081293400373,0.87081655452152,0.884256677500515,0.9005218061712,0.919464320861049,0.940463070117787,0.962165036846983,0.982170525852571,0.996664379730076,1,]
y = [0.496880137843736,0.892699355939761,1.31779053938845,1.76113078591939,2.2161404397406,2.68025205534181,3.15237712881927,3.63160088718474,4.11681195972416,4.60674953445045,5.10016093029361,5.59593557630099,6.09317746752204,6.59122126982431,7.08961104581689,7.58806045776256,8.08640803576273,8.58457519572663,9.08253011512446,9.58025775226833,10.0777349934011,10.5749095671634,11.0716819405241,11.5678904118526,12.0633008076957,12.5576032345323,13.0504183996379,13.5413128760906,14.0298073406874,14.5152905360297,15,]
from matplotlib import pyplot as plt
from numpy import *
x = matrix(x)
y = matrix(y)
xmin = -matrix(x)
y = matrix(y)
plt.plot(x.T, y.T,'-')
plt.plot(xmin.T, y.T,'-')
plt.axis('equal')
plt.show()
