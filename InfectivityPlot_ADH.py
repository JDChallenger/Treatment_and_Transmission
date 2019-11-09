import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import math
fig, ax = plt.subplots()

#this time make a vector for the time yourself!
L = 1
#time = np.arange(0, L, 0.025/24)
data = np.loadtxt('WHM_Results_gam_PK_ADH.txt')
dataP = np.loadtxt('WHM_Results_para_PK_ADH.txt')
dataPK = np.loadtxt('WHM_Results_AM_LMF_ADH.txt')
#print data[997]
#print len(time)
print len(dataP)

font4 = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

#gam = np.log10(data[1::3])
#gam2 = data[1::3]
#infect = data[2::3]
time = data[:,0]*(1.0/24.0)
print(time[9600-1])
gam = np.log10(data[:,1])
gam2 = data[:,1]
infect = data[:,2]

para = np.log10(dataP[4:204])
timeP = np.arange(0, 400, 2)
thresh = dataP[3]
start = dataP[0]*(0.2/24)
start2 = dataP[1]*(0.2/24)

if start2 == 0:
	inf_start = infect[round(0.2*dataP[0])] #Allow for delta changing with retreatment??
	inf_start2 = 0
else:
	inf_start = infect[round(0.2*dataP[0])]
	inf_start2 = infect[round(0.2*dataP[1])]
print(inf_start)
print(inf_start2)

timePK = (1.0/24.0)*dataPK[:,0]
am = dataPK[:,1]
dha = dataPK[:,2]
lmf = dataPK[:,3]
dlf = dataPK[:,4]

T = 20 #Default
for t in range(0, 9600):
    if gam2[t] > 0.00001 or para[math.floor(t/48)] > -5: #parasitaemia already 'logged'
    	T = time[t]#*(t);
print(T)
#T = 400

#mx = max(data[3:9603])
mxP = max(para)
listm = max(infect)
#mxl = max(listm)
#mxl = fn4(mx)

AUC = 0.0
for t in range(0, 9600-1):
	AUC += (1.0/24.0)*0.5*(infect[t] + infect[t+1])

print "AUC: ", AUC
AUC = round(AUC,2)

plt.figure(1)
plt.subplot(311)
plt.plot(time,gam, 'b-')#, t2, f(t2), 'k')
#plt.plot((1.0/24.0)*time,thresh, 'g-')
#plt.plot([0,T+3],[np.log10(thresh),np.log10(thresh)],'k--')
#plt.plot([0,T+3],[np.log10(50),np.log10(50)],'c-')
#plt.plot([(0.2/24.0)*data[0],(0.2/24.0)*data[0]],[-5.0,6],'k--')
plt.plot(timeP,para, 'r-')
plt.plot((start,start),(-5,mxP*1.2),'k--')
if start2 > 0 :
	plt.plot((start2,start2),(-5,mxP*1.2),'k--')
pl.xlim(0, T)
pl.ylim(-5.0, mxP*1.2)
plt.xticks([])
#pl.xlabel('Time (Days)')
#pl.ylabel(r'$\mathrm{log}_{10}[\mathrm{Parasites}/\mu L]$')
pl.ylabel('log P. Density')

plt.subplot(312)
#plt.plot((1.0/24.0)*time, fn(gam2), 'r--')
#plt.plot((1.0/24.0)*time, fn2(gam2), 'g-')
#plt.plot((1.0/24.0)*time, fn3(gam2), 'c-')
plt.plot(time, infect, 'k-')
plt.plot((start,start),(inf_start,1),'k--',zorder=1)
if start2 > 0 :
	plt.plot((start2,start2),(inf_start2,1),'k--',zorder=1)
#plt.plot([(0.2/24.0)*data[0],(0.2/24.0)*data[0]],[-5.0,6],'k--')
#plt.plot((1.0/24.0)*time,fn(gam2)
pl.xlim(0, T)
pl.ylim(0.0, 0.01+listm*1.15)
plt.xticks([])
pl.fill_between(time, 0, infect, facecolor='orange')
#pl.xlabel('Time (Days)')
pl.ylabel('P. (Infection)')
plt.text(0.85*T, 0.88*(0.008+listm*1.15), 'AUC: '+str(AUC),fontdict=font4)#bbox=dict(facecolor='red', alpha=0.5)
#plt.text(0.85*T, 0.72*(0.008+mxl*1.15), 'AUC2: '+str(AUC2),fontdict=font2)
#ax.annotate('local max', xy=(0.05, 0.1), xytext=(0.03, 0.05),
#           arrowprops=dict(facecolor='black', shrink=0.05),
#          )

#plt.subplot(413)
#plt.plot([(0.2/24.0)*data[0],(0.2/24.0)*data[0]],[0.0,1.05*max(am+dha)],'k--')
#plt.plot((1.0/24.0)*timePK, am, 'c-')
#plt.plot((1.0/24.0)*timePK, dha, 'k-')
#plt.locator_params(numticks=5)
#ax.set_xticks(ax.get_xticks()[::2])
#pl.xlabel('Time (Days)')
#plt.xticks([])
#pl.ylabel('AM/DHA Conc.')
#plt.plot([(0.2/24.0)*data[0],(0.2/24.0)*data[0]],[0.0,1.05*max(am+dha)],'k--')
#pl.xlim(0, T)#min([T+3,(0.2/24.0)*data[0]+19]))

plt.subplot(313)
plt.plot(timePK, lmf, 'g-')
plt.plot((start,start),(0,1.05*max(lmf)),'k--')
if start2 > 0 :
	plt.plot((start2,start2),(0,1.05*max(lmf)),'k--')
#plt.plot((1.0/24.0)*timePK, dlf, 'r-')
#plt.plot([(0.2/24.0)*data[0],(0.2/24.0)*data[0]],[0.0,1.05*max(lmf)],'k--')
pl.xlabel('Time (Days)')
pl.ylabel('LMF Conc.')
pl.ylim(0.0, 1.05*max(lmf))
pl.xlim(0, T)#min([T+3,(0.2/24.0)*data[0]+19]))
ax.set_yticks(ax.get_yticks()[::2])

plt.savefig("GamDensityAndInfectivityADH.pdf", format='pdf')
plt.show()