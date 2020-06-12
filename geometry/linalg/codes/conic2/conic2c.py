import numpy as np
import matplotlib.pyplot as plt
import subprocess
import  shlex
import math

#coeffecients
a=2
b=7
c=3
#writing the equation
x = np.linspace(-4,4,100)
p=np.poly1d([a,b,c])
for i in range(len(x)):
    y=p(x)

#zeroes
s=math.sqrt(b**2-4*a*c)
k1=(-1*b+s)/(2*a)
k2=(-1*b-s)/(2*a)


A=np.array([k1,p(k1)])
B=np.array([k2,p(k2)])
plt.plot(A[0], A[1], 'o')
plt.text(A[0] * (1 + 0.001), A[1] * (1 - 0.01) , 'k1')
plt.plot(B[0], B[1], 'o')
plt.text(B[0] * (1-0.005), B[1] * (1-0.01) , 'k2')
plt.plot(x,y,label='$p(x)=2x^2+7x+3$')
plt.xlabel('$x$')
plt.ylabel('$p(x)$')
plt.legend(loc='best')
plt.grid() # minor
#plt.axis('equal')
plt.savefig('./pyfigs/conic2c.eps')
plt.show() 
