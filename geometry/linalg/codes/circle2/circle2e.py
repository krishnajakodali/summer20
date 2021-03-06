import numpy as np
import matplotlib.pyplot as plt
from coeffs import *
import subprocess
import shlex
import math
#parameters
a=5
b=4
#centre and radius of the circle
O=np.array([-a,-b])
r=math.sqrt(a**2-b**2)
#point on the circle
A=np.array([O[0]+r,O[1]])

#Generating all lines
x_OA = line_gen(O,A)

#Plotting all lines
plt.plot(x_OA[0,:],x_OA[1,:],label='$r=3$')
#plt.plot(x_OB[0,:],x_OB[1,:])
#plotting circle
len=100
theta = np.linspace(0,2*np.pi,len)
x_circ = np.zeros((2,len))
x_circ[0,:] = r*np.cos(theta)
x_circ[1,:] = r*np.sin(theta)
x_circ = (x_circ.T + O).T
plt.plot(x_circ[0,:],x_circ[1,:],label='$circle(e)$')


#plotting points
plt.plot(O[0], O[1], 'o')
plt.text(O[0] * (1 + 0.001), O[1] * (1 - 0.001) , 'O')

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='upper right')
plt.grid() # minor
plt.axis('equal')
plt.savefig('./pyfigs/circle2e.eps')
plt.show()

