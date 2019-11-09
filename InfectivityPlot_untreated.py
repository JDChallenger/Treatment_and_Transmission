import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import math
fig, ax = plt.subplots()

font4 = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

data=np.loadtxt('WHM_Gametocytes_Results.txt')
b=0
for t in range(0, 399):
    if data[t,2] > 0.00001 or data[t,1] > 0.00001:
    	T = t;
    if data[t,1] > b:
    	b = data[t,1]
    if data[t,1] < 0.00001 and data[t,2] < 0.00001:
    	break

print(b)
print(T)

mxl = max(data[:,3])

AUC = 0
for t in range(0, 398):
	AUC += 0.5*2*(data[t,3]+data[t+1,3]) #time step is 2 days
#Round to something sensible	
AUC =round(AUC,2)	

plt.figure(1)
plt.subplot(211)
pl.plot(data[:,0],data[:,1])#It seems we can add as many plots here as we want!
pl.plot(data[:,0],data[:,2])#Gametocytes
pl.plot([0,(2 * T) + 6],[40,40])#Adds the microscopy threshold. What value?
pl.xlabel('Time (Days)')
pl.ylabel('Parasitaemia')
pl.yscale('log')
plt.xticks([])
pl.xlim(0.0, 6 + (2 * T)) #Increased now, to allow for mature gametocytes to circulate
pl.ylim(0.0001, b * 1.25) #Scaling of 1.25 just provides some space above 1st peak

plt.subplot(212)
pl.plot(data[:,0],data[:,3],'k-')
pl.xlabel('Time (Days)')
pl.ylabel('Probability of Transmission')#?
pl.fill_between(data[:,0],0, data[:,3], facecolor='orange')
plt.text(0.81*2*T, 0.85*(0.01+mxl*1.05), 'AUC: '+str(AUC),fontdict=font4)
pl.xlim(0.0, 6 + (2 * T)) #Increased now, to allow for mature gametocytes to circulate
pl.ylim(0.00001, 0.01+mxl*1.15) 
ax.set_yticks(ax.get_yticks()[::2])

plt.savefig("WHM_Gametocytes_one_run.pdf", format='pdf')
plt.show()