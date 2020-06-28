import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import bernoulli
simlen=int(1e5)
sample_size=5
#queen from 5 cards
event_size1=1
prob1=event_size1/sample_size
#ace from 4 cards
event_size2=1
prob2=event_size2/(sample_size-1)
#queen from 4 cards
event_size3=0
prob3=event_size3/(sample_size-1)
#Simulations using Bernoulli r.v.
data_bern1 = bernoulli.rvs(size=simlen,p=prob1)
data_bern2= bernoulli.rvs(size=simlen,p=prob2)
data_bern3 = bernoulli.rvs(size=simlen,p=prob3)



#Calculating the number of favourable outcomes
err_ind1 = np.nonzero(data_bern1 == 1)
err_ind2 = np.nonzero(data_bern2 == 1)
err_ind3 = np.nonzero(data_bern3 == 1)
print(err_ind1)
print(err_ind2)
print(err_ind3)

#calculating the probability
err_n1 = np.size(err_ind1)/simlen
err_n2 = np.size(err_ind2)/simlen
err_n3 = np.size(err_ind3)/simlen


#Theory vs simulation
print(err_n1,prob1)
print(err_n2,prob2)
print(err_n3,prob3)
