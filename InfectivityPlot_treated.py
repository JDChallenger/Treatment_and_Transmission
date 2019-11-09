import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import math
fig, ax = plt.subplots()

#this time make a vector for the time yourself!
L = 1
#time = np.arange(0, L, 0.025/24)
data = np.loadtxt('WHM_Results_gam_and_infect_PK_ind.txt')
dataP = np.loadtxt('WHM_Results_para_PK_ind.txt')
#print data[997]
#print len(time)
#print len(dataP)

font4 = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

T = 20 #Default
for t in range(0, 9600):
    if data[t,0] > 0.00001:
    	T = (1.0/24.0)*(t-2);
#print(T)


gam = np.log10(data[:,0])
infect = data[:,1]
#gam2 = data[3:9603]
time = np.arange(0, 9600, 1)

para = np.log10(dataP[4:9604])
timeP = np.arange(0, 400, 2)
thresh = dataP[3]
start = dataP[0]*(0.2/24)
start2 = dataP[1]*(0.2/24)
#print(dataP[0])
#print(dataP[1])
#print(start)
#print(start2)
if start2 == 0:
	inf_start = infect[round(0.2*dataP[0])] #Allow for delta changing with retreatment??
	inf_start2 = 0
else:
	inf_start = infect[round(0.2*dataP[0])]
	inf_start2 = infect[round(0.2*dataP[1])]
print(inf_start)
print(inf_start2)

mx = max(data[:,0])
mxP = max(para)
mxl = max(infect)
#mxl = max(listm)

#print len(timeP)
#print len(para)
#print thresh

AUC = 0.0
for t in range(0, 9600-1):
	AUC += (1.0/24.0)*0.5*(infect[t]+infect[t+1])
print "AUC: ", AUC

AUC =round(AUC,2)

plt.figure(1)
plt.subplot(211)
plt.plot(timeP,para, 'r-',label="Asexual Parasites")
plt.plot((1.0/24.0)*time,gam, 'b-',label="Gametocytes")#, t2, f(t2), 'k')
#plt.plot((1.0/24.0)*time,thresh, 'g-')
#plt.plot((0,T+3),(np.log10(thresh),np.log10(thresh)))
plt.plot((start,start),(-5,mxP*1.2),'k--')
if start2 > 0 :
	plt.plot((start2,start2),(-5,mxP*1.2),'k--')
#plt.plot([0,T+3],[np.log10(thresh),np.log10(thresh)],'k--')
#plt.plot([(0.2/24.0)*dataP[0],(0.2/24.0)*data[0]],[-5.0,6],'k--')
pl.xlim(0, T+3)
pl.ylim(-5.0, mxP*1.2)
pl.ylabel('log10[Parasite Density]')
plt.legend(frameon=False)

plt.subplot(212)
#plt.plot((1.0/24.0)*time, fn(gam2), 'r--')
plt.plot((1.0/24.0)*time, infect, 'k-',zorder=10)
pl.xlim(0, T+3)
pl.ylim(0.0, 0.01+mxl*1.15)
pl.xlabel('Time (Days)')
pl.ylabel('Probability of Infection')
plt.plot((start,start),(inf_start,1),'k--',zorder=1)
if start2 > 0 :
	plt.plot((start2,start2),(inf_start2,1),'k--',zorder=1)
pl.fill_between((1.0/24.0)*time,0, infect, facecolor='orange', alpha=1.0,zorder=5)
plt.text(0.85*T, 0.85*(0.01+mxl*1.15), 'AUC: '+str(AUC),fontdict=font4)#bbox=dict(facecolor='red', alpha=0.5)
#ax.annotate('local max', xy=(0.05, 0.1), xytext=(0.03, 0.05),
 #           arrowprops=dict(facecolor='black', shrink=0.05),
  #          )
plt.savefig("GamDensityAndInfectivityRetreat.pdf", format='pdf')
plt.show()