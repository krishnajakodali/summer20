import numpy as np
import matplotlib.pyplot as plt
import math
import cmath
import random
#if using termux
import subprocess
import shlex
from scipy import signal


TimeSlot=2e-3 #Transmit time duration
SNR = 18 #Signal to noise ratio
Rs = 185e3 # symbol rate
a=np.array([1+0j,1/math.sqrt(2)+1j*1/math.sqrt(2) ,1j ,-1/math.sqrt(2)+1j*1/math.sqrt(2), -1 ,-1/math.sqrt(2)-1j*1/math.sqrt(2), -1j ,1/math.sqrt(2)-1j*1/math.sqrt(2)])

Ak=[]
for i in range(10000):
	Ak.append(0)

for ii in range(10000):
	Ak[ii] = a[random.randint(0,7)]
ar=[x.real for x in a]
ai=[x.imag for x in a]
plt.scatter(ar,ai)
plt.title('8-PSK')
plt.grid()
plt.show()


# Channel creation and channel modelling
Rsym = Rs; M = 8;                  # Input symbol rate
Rbit = Rsym * math.log2(M);      #Input bit rate
Nos = 1;                    # Oversampling factor
ts = (1/Rsym) / Nos;
pg=np.array([0.9417 - 0.2214j ,-0.1596 - 0.0874j ,-0.0644 - 0.0163j, -0.0645 - 0.0387j, -0.0751 + 0.0467j])

#readsfrom tex file command File_object.readline(n)
#pg=dlmread('path_gains.dat',',',[0,0,0,4])
pd=np.array([0, 2.0000e-06 ,4.0000e-06, 6.0000e-06, 8.0000e-06])/ts
g=[]
for n in range(0,1600):
     g.append(0)
     for k in range(5):
         g[n]=g[n]+pg[k]*np.sinc(pd[k]-n+800)
Rk=np.convolve(Ak,g,'same')
noise = (1/math.sqrt(2))*(np.random.randn(len(Rk)) + 1j*np.random.randn(len(Rk))) #Initial noise vector
P_s =np.var(Rk)  #Signal power
P_n = np.var(noise)  #Noise power
# Defining noise scaling factor based on the desired SNR:
noise_scaling_factor = math.sqrt((P_s/P_n)/10**(SNR/10)) 
Rk_noisy=Rk+noise*noise_scaling_factor # Received signal
print(len(Rk_noisy))
rkr=[x.real for x in Rk]
rki=[x.imag for x in Rk]
plt.scatter(rkr,rki)
plt.title('Recieved constellation')
plt.grid()
plt.show()

Ek=[]
hTap=11 #Channel Taps
beta = 0.001 # step-size of the algorithm
c_LMS = np.zeros([hTap,1]) #equalizer coefficients, initializations

for i in range(0,len(Rk_noisy)-int((hTap-1)/2)-1):
	Ek.append(0)
for i in range(int((hTap-1)/2),len(Rk_noisy)-int((hTap-1)/2)-1):
	rkt =Rk_noisy[i-int((hTap-1)/2):i+int((hTap-1)/2)+1] 
	rk=np.flipud(rkt.reshape(-1,1)) # recieved signal vector
	Ek[i] = Ak[i] - np.matmul(np.transpose(c_LMS),rk) #Error signal, we assume a known symbol sequence
	c_LMS = np.add(c_LMS,beta*Ek[i]*(np.conj(rk))) # LMS update !

# MSE Equalizer%%%%%%%%%%%%%%%%%%%%%%%LMS Algorithm%%%%%%%%%%%%%%%%
# Initialization
FII = np.zeros([hTap,hTap]) # Autocorrelation Matrix initialization
alfa = np.zeros([hTap,1]) # Cross-correlation vector initialization
# Estimating FII and alfa using sample estimates based on training symbols
# Notice that here we use all the generated data as training symbols
for i in range(int((hTap-1)/2) , len(Rk_noisy)-int((hTap-1)/2)-1): #16:length(Rk_noisy)-15, 
	rkt =Rk_noisy[i-int((hTap-1)/2) :i+int((hTap-1)/2)+1] # Received signal vector
	rk=np.flipud(rkt.reshape(-1,1))
	FII =np.transpose(FII + rk*np.conj(rk)) # Autocorrelation matrix
	alfa = alfa + (Ak[i]*(np.conj(rk)))
# FII
FII = FII/(len(Rk_noisy)-(hTap-1)) # Final sample estimate of the autocorrelation matrix
alfa = alfa/(len(Rk_noisy)-(hTap-1)) # Final sample estimate of the cross-correlation vector
c_MSE = np.matmul(np.linalg.inv(np.conj(FII)),alfa) #Equalizer coefficients
#
# Plotting
y_LMS=signal.lfilter(np.ndarray.flatten(c_LMS),1,np.ndarray.flatten(Rk))
y_MSE=signal.lfilter(np.ndarray.flatten(c_MSE),1,np.ndarray.flatten(Rk))

Lr=[x.real for x in y_LMS]
Li=[x.imag for x in y_LMS]

Mr=[x.real for x in y_MSE]
Mi=[x.imag for x in y_MSE]


plt.scatter(Lr,Li)
plt.title('LMS equallized constellation')
plt.grid()
plt.show()

plt.scatter(Mr,Mi)
plt.title('MSE equalized constellation')
plt.grid()
plt.show()
