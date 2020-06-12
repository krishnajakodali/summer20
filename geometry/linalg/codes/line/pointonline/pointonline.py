
import numpy as np
import matplotlib.pyplot as plt
from coeffs import *
#using termux
import subprocess
import  shlex

x=np.linspace(-5,5,100)
y1=(12-4*x)/3
y2=(-2*x)/5
y3=4/3+x*0

plt.plot(x,y1,label='4x+3y=12')
plt.plot(x,y2,label='2x+5y=0')
plt.plot(x,y3,label='3y=4')
plt.title('Lines representing the three equations')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor


#different points on a line
A1 = np.array([0,4])
B1= np.array([1,8/3])
A2 = np.array([0,0])
B2 = np.array([1,-2/5])
A3 = np.array([0,4/3])
B3 = np.array([1,4/3])

#plotting point
plt.plot(A1[0], A1[1], 'o')
plt.text(A1[0] * (1 - 0.1), A1[1] * (1 + 0.1) , 'A1')
plt.plot(B1[0], B1[1], 'o')
plt.text(B1[0] * (1 - 0.1), B1[1] * (1 + 0.1) , 'B1')
plt.plot(A2[0], A2[1], 'o')
plt.text(A2[0] * (1 - 0.1), A2[1] * (1 + 0.1) , 'A2')
plt.plot(B2[0], B3[1], 'o')
plt.text(B2[0] * (1 - 0.1), B2[1] * (1 + 0.1) , 'B2')
plt.plot(A3[0], A3[1], 'o')
plt.text(A3[0] * (1 - 0.1), A3[1] * (1 + 0.1) , 'A3')
plt.plot(B3[0], B3[1], 'o')
plt.text(B3[0] * (1 - 0.1), B3[1] * (1 + 0.1) , 'B3')

plt.savefig('./pyfigs/pointonline.eps')

#print(3)
plt.show()
#print(4)




