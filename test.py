import numpy as np
import matplotlib.pyplot as plt

time = np.linspace(0, 10, 100)
sinout = np.sin(time)
cosout = np.cos(time)

plt.plot(time, sinout, label='sin(t)')
plt.plot(time, cosout, label='cos(t)')
plt.grid()
plt.legend()
plt.show()

