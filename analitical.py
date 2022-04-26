import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

L=100
h=1
a=5e0
q=(a*np.pi/L)**2

X = np.arange(0, L+h, h)
T = np.arange(0, L+h, h)
Z = np.zeros(int(L/h)+1)
X, T = np.meshgrid(X, T)

for i in range(1, 1000+1, 2):
    z_t=np.exp(-q*T*i*i)
    z_x=np.sin(np.pi*i*X/L)
    Z=Z+z_t*z_x*400/(np.pi*i)

fig=plt.figure()
ax=fig.add_subplot(projection='3d')
surf = ax.plot_surface(X, T, Z, cmap=cm.coolwarm, alpha = 0.5,
                       linewidth=0, antialiased=False)


plt.show()
