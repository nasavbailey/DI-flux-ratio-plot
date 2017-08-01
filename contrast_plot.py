#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 11:33:07 2017

@author: meshkat
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})
matplotlib.rcParams.update({'font.family':'Times New Roman'})
matplotlib.rcParams.update({'font.size': 10})


fig=plt.figure()
ax1=fig.add_subplot(111)

###Define path where to find data and where to save to
path='/Users/meshkat/python/WFIRST/'


#########################################################################
#### ELT Section
#### Note, comment out this whole section if you don't want to include the ELTs.
#range_x=np.array((0.03,1))
#pessimistic_y=np.array((10**-5,10**-9))
#optimistic_y=np.array((10**-8,10**-9))
#ax1.plot(range_x,pessimistic_y,color='pink', linestyle='--',linewidth=1)
#ax1.plot(range_x,optimistic_y,color='pink', linestyle='--',linewidth=1)
#plt.fill_between(range_x, pessimistic_y, optimistic_y, color='pink', alpha='0.2')
#ax1.plot(np.array((0.03,1)),np.array((10**-6,10**-9)),color='red', linestyle='--',linewidth=1)
#
####Now add text
#plt.text(0.1,2*10**-7,'Future ELTs?',color='Red',horizontalalignment='left',verticalalignment='top',rotation=-19)

#########################################################################
###Planet contrasts at certain separations ((sep in arcsec, contrast))
HR8799b=np.array((1.4,8.5E-6))
HR8799c=np.array((0.9,10**-5))
HR8799d=np.array((0.6, 2.9E-5))
HR8799e=np.array((0.35,3E-5))
betaPicb=np.array((0.4,10**-4))
HD95086b=np.array((0.6,7E-6))
Eri51b=np.array((0.45,1E-6))
ProximaCenb=np.array((0.035,4E-8))

####Plot self luminous planets in orange


markersize_self_luminous=4
color_self_luminous='red'

plt.plot(HR8799b[0],HR8799b[1],marker='o',color=color_self_luminous,markersize=markersize_self_luminous)
plt.plot(HR8799c[0],HR8799c[1],marker='o',color=color_self_luminous,markersize=markersize_self_luminous)
plt.plot(HR8799d[0],HR8799d[1],marker='o',color=color_self_luminous,markersize=markersize_self_luminous)
plt.plot(HR8799e[0],HR8799e[1],marker='o',color=color_self_luminous,markersize=markersize_self_luminous)
plt.plot(betaPicb[0],betaPicb[1],marker='o',color=color_self_luminous,markersize=markersize_self_luminous)
plt.plot(HD95086b[0],HD95086b[1],marker='o',color=color_self_luminous,markersize=markersize_self_luminous)
plt.plot(Eri51b[0],Eri51b[1],marker='o',color=color_self_luminous,markersize=markersize_self_luminous)
plt.plot(ProximaCenb[0],ProximaCenb[1],marker='o',color=color_self_luminous,markersize=markersize_self_luminous)


###Plot names of planets
plt.text(HR8799b[0],HR8799b[1],' HR 8799 b ',color='black',horizontalalignment='left',verticalalignment='center')
plt.text(HR8799c[0],HR8799c[1],' HR 8799 c ',color='black',horizontalalignment='right',verticalalignment='center')
plt.text(HR8799d[0],HR8799d[1],' HR 8799 d ',color='black',horizontalalignment='left',verticalalignment='center')
plt.text(HR8799e[0],HR8799e[1],' HR 8799 e ',color='black',horizontalalignment='right',verticalalignment='center')
plt.text(betaPicb[0],betaPicb[1],' beta Pic b ',color='black',horizontalalignment='left',verticalalignment='center')
plt.text(HD95086b[0],HD95086b[1],' HD 95086 b ',color='black',horizontalalignment='right',verticalalignment='top')
plt.text(Eri51b[0],Eri51b[1],' 51 Eri b ',color='black',horizontalalignment='left',verticalalignment='center')
plt.text(ProximaCenb[0],ProximaCenb[1],'  Proxima Cen b ',color='black',horizontalalignment='left',verticalalignment='center')

#########################################################################
#########Read in contrast curves for Exoplanet Mode and Disk Mode    
a_e = np.loadtxt(path+'exoplanet_mode.csv',delimiter=',')
a_d = np.loadtxt(path+'disk_mode.csv',delimiter=',')


arcsec_exoplanet=a_e[:,1]
contrast_exoplanet=a_e[:,2]
arcsec_disk=a_d[:,1]
contrast_disk=a_d[:,2]

plt.plot(arcsec_exoplanet,contrast_exoplanet,color='black',linewidth=2)
plt.plot(arcsec_disk,contrast_disk,color='black',linewidth=2)

#########################################################################
########Add text stating where the CGI Baseline science requirements and CGI threshold technical requirements are
plt.text(0.15,9*10**-10,'CGI Baseline Science Requirement',color='black',horizontalalignment='left',verticalalignment='top')


######Add Technical requirement line and text
#plt.plot(np.array((0.15,1.)),np.array((10**-7,10**-7)),color='black',linewidth=2)
#plt.text(0.15,5*10**-8,'CGI Threshold Technical Requirement',color='black',horizontalalignment='left',verticalalignment='top')

#########################################################################
########Add HST, JWST, SPHERE, GPI contrast curves
a_HST = np.loadtxt(path+'hst_1.txt')
a_JWST = np.loadtxt(path+'jwst_1.txt')
a_SPHERE = np.loadtxt(path+'SPHERE_Vigan.txt')
a_GPI = np.loadtxt(path+'GPIES_T-type_contrast_curve_2per.txt')

##Read the specific arcsec and contrast columns we need
arcsec_HST=a_HST[:,1]
contrast_HST=a_HST[:,2]
arcsec_JWST=a_JWST[:,1]
contrast_JWST=a_JWST[:,2]
arcsec_SPHERE=a_SPHERE[:,0]
contrast_SPHERE=10**(-a_SPHERE[:,1]/2.5) ###Note this is in magnitudes, need to convert to flux
arcsec_GPI=a_GPI[:,0]
contrast_GPI=a_GPI[:,1]

###Plot the contrast curves
plt.plot(arcsec_HST,contrast_HST,color='black',linewidth=1)
plt.plot(arcsec_JWST,contrast_JWST,color='red',linewidth=1)
plt.plot(arcsec_SPHERE,contrast_SPHERE,color='red',linewidth=1)
plt.plot(arcsec_GPI,contrast_GPI,color='red',linewidth=1)

###Label the curves
plt.text(0.2,5*10**-7,'SPHERE-Sirius',color='red',horizontalalignment='left',rotation=-20,fontsize=10)
plt.text(0.19,2.5*10**-6,'GPI',color='red',horizontalalignment='left',rotation=-20,fontsize=10)
plt.text(1.4,7*10**-7,'JWST NIRCam',color='red',horizontalalignment='left',rotation=-30,fontsize=10)
plt.text(2.2,10**-8,'HST ACS',color='black',horizontalalignment='left',rotation=-30,fontsize=10)


#########################################################################
###Add Earth and Jupiter
Earth=np.array((0.09,1*10**-10))
Jupiter=np.array((0.51, 10**-9))
plt.plot(Earth[0],Earth[1],marker='o',color='blue',markersize=markersize_self_luminous)
plt.plot(Jupiter[0],Jupiter[1],marker='o',color='#FF7800',markersize=markersize_self_luminous)
plt.text(Earth[0],Earth[1],' Earth ',color='blue',horizontalalignment='left',verticalalignment='center')
plt.text(Jupiter[0],Jupiter[1],' Jupiter ',color='#FF7800',horizontalalignment='left',verticalalignment='center')

#########################################################################
###Add RV planets
a_RV = np.loadtxt(path+'RVtable.txt')
arcsec_RV=a_RV[:,0]
contrast_RV=a_RV[:,1]
plt.plot(arcsec_RV,contrast_RV, 'ko', markersize=markersize_self_luminous)




#########################################################################
###Plot axes, tick mark ajdusting, legend, etc.
x_ticklabels=np.array((0.1,1))##arcsec  
x_ticks=x_ticklabels
x_ticks_minor=np.array((0.05,0.5))
x_ticklabels_minor=x_ticks_minor

y_ticks=np.array((10**-11,10**-10,10**-9,10**-8,10**-7,10**-6,10**-5,10**-4,10**-3))
y_ticklabels=np.log10(y_ticks)
##


ax1.set_ylabel('Flux raio to host star')
ax1.set_xlabel('Separation [arcsec]') 

ax1.set_ylim(10**-11,10**-2.9)
ax1.set_xlim(0.03,4)
ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.set_xticks(x_ticks)
ax1.set_xticklabels(x_ticklabels)


###Add grid lines
plt.grid(b=True, which='major', color='b', linestyle='--', alpha=0.1)

####Add custom legend (probably a better way to do this?)
plt.plot(np.array((0.035,0.15)),np.array((4*10**-4,4*10**-4)), 'k--',linewidth=1, alpha=0.5)
plt.plot(np.array((0.035,0.15)),np.array((4*10**-7,4*10**-7)), 'k--',linewidth=1, alpha=0.5)
plt.plot(np.array((0.035,0.035)),np.array((4*10**-7,4*10**-4)), 'k--',linewidth=1, alpha=0.5)
plt.plot(np.array((0.15,0.15)),np.array((4*10**-7,4*10**-4)), 'k--',linewidth=1, alpha=0.5)
plt.plot(0.039,2*10**-4,marker='o',color='black',markersize=markersize_self_luminous)
plt.text(0.039,2*10**-4,'  Known RV (visible)',color='black',horizontalalignment='left',verticalalignment='center',fontsize=8)
plt.plot(0.039,9*10**-5,marker='o',color=color_self_luminous,markersize=markersize_self_luminous)
plt.text(0.039,9*10**-5,'  Known self-luminous (NIR)',color=color_self_luminous,horizontalalignment='left',verticalalignment='center',fontsize=8)
plt.text(0.039,10**-5,'Reflected light contrast \nas seen at quadrature',color='black',horizontalalignment='left',verticalalignment='center',fontsize=10)
plt.text(0.039,1*10**-6,'Solar system planets \nas seen from 10 pc',color='black',horizontalalignment='left',verticalalignment='center',fontsize=10)

###Add minor tick labels
#ax1.set_yticks(y_ticks)
#ax1.set_yticklabels(y_ticklabels)

#ax1.set_xticks(x_ticks_minor, minor = True)
#ax1.set_xticklabels(x_ticklabels_minor, minor=True)

plt.savefig(path+'flut_ratio_plot.pdf')
