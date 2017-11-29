#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: meshkat (original; July 31, 2017)
November 2017: VPB reorganized & updated contrast curves
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import os
from matplotlib import rcParams
from astropy.io import ascii
#import decimal
from matplotlib.pyplot import gca

### User-defined options

# Which contrast curves to include?
include_ELT     = False
include_HABEX   = False
include_ACS     = False
include_NICMOS  = True
include_STIS    = True
include_NIRCAM  = True
include_ACS     = False # don't use until verified
include_SPHERE  = True
include_GPI     = True
include_DI_H    = True # real H-band contrasts of known directly imaged planets
include_DI_750_extrap = True # COND/BT-Settl model extrapolations to ~750nm
include_DI_550_extrap = True # COND/BT-Settl model extrapolations to ~550nm

color_by_lambda = True # colorcode contrast curve lines by wavelength?

###Define path where to find data and where to save to
path = '' # leave blank if this script is in the same folder as the data (default)


########################################################################
### Plot setup

rcParams.update({'figure.autolayout': True})
matplotlib.rcParams.update({'font.family':'Times New Roman'})
matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams['mathtext.fontset'] = 'stix'

fig=plt.figure()
ax1=fig.add_subplot(111)
#ax2 = ax1.twinxy()

markersize_points=4
ccfs = 8 # contrast curve font size
ccc = 'darkviolet' # default contrast curve color
cclw = 2 # default contrast curve line width

if color_by_lambda:
    c_550 = 'cadetblue'
    c_750 = 'goldenrod'
    c_yjh = 'coral'
    c_k = 'firebrick'
    c_h   = 'red'
    c_l   = 'darkred'

    plt.plot([1,1],[1,1],color=c_550,linewidth=cclw, label='Visible')
    plt.plot([1,1],[1,1],color=c_yjh,linewidth=cclw, label='YJH-band')
    plt.plot([1,1],[1,1],color=c_h,linewidth=cclw, label='H-band')
    plt.plot([1,1],[1,1],color=c_k,linewidth=cclw, label='K-band')

else:
    c_550 = ccc
    c_750 = ccc
    c_yjh = ccc
    c_k = ccc
    c_h   = ccc
    c_l   = ccc

########################################################################
# auto-generated caption. See README for how to comment datafiles.

caption = '** Please see individual datafiles for full descriptions.** \n\n'

def extract_short_caption(filename):
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    for l in lines:
        if '#short caption:' in l.lower():
            return l.split('caption:')[1].strip()
    # if no caption in text file
    print '\n**** WARNING **** no caption for '+filename+'\n'
    return ''


########################################################################
# optional suffix, depending on what's plotted?
file_name_end = ''

########################################################################
### ELT

if include_ELT:
    #file_name_end += '_ELT'
    range_x=np.array((0.03,1))
    pessimistic_y=np.array((10**-5,10**-8))
    optimistic_y=np.array((10**-8,10**-9))
    ax1.plot(range_x,pessimistic_y,color='pink', linestyle='--',linewidth=cclw)
    ax1.plot(range_x,optimistic_y,color='pink', linestyle='--',linewidth=1)
    plt.fill_between(range_x, pessimistic_y, optimistic_y, color='pink', alpha='0.2')
    #ax1.plot(np.array((0.03,1)),np.array((10**-6,10**-9)),color='red', linestyle='--',linewidth=1)

    ###Now add text
    plt.text(0.1,2.5*10**-7,'Future ELTs (NIR)?',color='Red',horizontalalignment='left',\
    	verticalalignment='top',rotation=-19,fontsize=ccfs)


#########################################################################
### HabEx "goal" contrast curve

if include_HABEX:
    ax1.plot(np.array((.05,1.7)),np.array((10**-10,10**-10)),color='black', linestyle='--',linewidth=1)
    plt.text(0.2,9*10**-11,'HabEx Goal',color='black',horizontalalignment='left',\
    	verticalalignment='top',fontsize=ccfs)


#########################################################################
### NIRCAM contrast curve

if include_NIRCAM:
    fname = path+'jwst_nircam.txt'
    a_JWST = ascii.read(fname)
    plt.plot(a_JWST['arcsec'],a_JWST['210_contr'],color=c_k,linewidth=cclw-0.5,linestyle='--',label='')
    plt.text(0.9,2*10**-6,'JWST NIRCam',color=c_k,horizontalalignment='left',rotation=-35,fontsize=ccfs)
    caption += extract_short_caption(fname)+'\n'

#########################################################################
### NICMOS contrast curve
if include_NICMOS:
    fname = path+'HST_NICMOS_Min.txt' #path+'HST_NICMOS_Median.txt'
    a_NICMOS = ascii.read(fname)
    plt.plot(a_NICMOS['Rho(as)'],a_NICMOS['F160W_contr'],color=c_h,\
        linewidth=cclw,linestyle='-.',label='')
    plt.text(1.7,1*10**-6,'HST NICMOS',color=c_h,horizontalalignment='left',rotation=-20,fontsize=ccfs)
    caption += extract_short_caption(fname)+'\n'

#########################################################################
### STIS Bar5 contrast curve
if include_STIS:
    fname = path+'HST_STIS.txt'
    a_STIS = ascii.read(fname)
    plt.plot(a_STIS['Rho(as)'],a_STIS['KLIP_Contr'],color=c_550,\
        linewidth=cclw,linestyle='-.',label='')
    plt.text(0.2,5*10**-5,'HST STIS',color=c_550,horizontalalignment='left',rotation=-40,fontsize=ccfs)
    caption += extract_short_caption(fname)+'\n'

#########################################################################
### ACS contrast curve

if include_ACS:
	print "ACS data is not verified. Skipping"
	#a_HST = np.loadtxt(path+'hst_1.txt')
	#arcsec_HST=a_HST[:,1]
	#contrast_HST=a_HST[:,2]
	#plt.plot(arcsec_HST,contrast_HST,color='black',linewidth=1)
	#plt.text(2.2,10**-8,'HST ACS',color='black',horizontalalignment='left',rotation=-30,fontsize=ccfs)


#########################################################################
### SPHERE contrast curve
if include_SPHERE:
    fname = path+'SPHERE_Vigan.txt'
    a_SPHERE = ascii.read(fname)
    a_SPHERE['Contrast'] = 10**(-0.4*a_SPHERE['delta'])

    # manually split into IFS and IRDIS, at 0.7", as per documentation.
    idx_yjh = a_SPHERE['Rho(as)'] <= 0.7 # IFS YJH
    idx_k12 = a_SPHERE['Rho(as)'] >= 0.7  # IRDIS K1-K2
    plt.plot(a_SPHERE['Rho(as)'][idx_yjh], a_SPHERE['Contrast'][idx_yjh], \
        color=c_yjh, linewidth=cclw, label='')
    plt.plot(a_SPHERE['Rho(as)'][idx_k12], a_SPHERE['Contrast'][idx_k12], \
        color=c_k, linewidth=cclw, label='')
    plt.text(0.15,1*10**-6,'SPHERE IFS',color=c_yjh,horizontalalignment='left',rotation=-20,fontsize=ccfs)
    plt.text(1.5,8*10**-8,'SPHERE IRDIS',color=c_k,horizontalalignment='left',rotation=-25,fontsize=ccfs)
    caption += extract_short_caption(fname)+'\n'

#########################################################################
### GPI H-band
if include_GPI:
    fname = path+'GPIES_T-type_contrast_curve_2per.txt'
    a_GPI = ascii.read(fname)
    plt.plot(a_GPI['Rho(as)'],a_GPI['H_contr'],color=c_h,linewidth=cclw,label='')
    plt.text(0.17,1*10**-5,'GPI',color=c_h,horizontalalignment='left',rotation=-20,fontsize=ccfs)
    caption += extract_short_caption(fname)+'\n'

#########################################################################
### WFIRST
a_e = np.loadtxt(path+'exoplanet_mode.csv',delimiter=',')
a_d = np.loadtxt(path+'disk_mode.csv',delimiter=',')


arcsec_exoplanet=a_e[:,1]
contrast_exoplanet=a_e[:,2]
arcsec_disk=a_d[:,1]
contrast_disk=a_d[:,2]

plt.plot(arcsec_exoplanet,contrast_exoplanet,color='black',linewidth=2)
plt.plot(arcsec_disk,contrast_disk,color='gray',linewidth=2) ###This is commented for now, can be added later

#########################################################################
########Add text stating where the CGI Baseline science requirements and CGI threshold technical requirements are
#plt.text(0.1,8*10**-10,'WFIRST CGI Baseline Exoplanet Detectability',color='black',horizontalalignment='left',verticalalignment='top',fontsize=10)


######Add Technical requirement line and text
#plt.plot(np.array((0.15,1.)),np.array((10**-7,10**-7)),color='black',linewidth=2)
#plt.text(0.15,5*10**-8,'CGI Threshold Technical Requirement',color='black',horizontalalignment='left',verticalalignment='top')



#########################################################################
###### -------------------- planets -------------------------  ##########
#########################################################################

#########################################################################
### Self luminous directly imaged planets

if include_DI_H or include_DI_extrap:
    a_DI = ascii.read(path+'DIplanets.txt')
    a_DI['547m_contr'] = 10**(a_DI['547m_delta']/-2.5)
    a_DI['763m_contr'] = 10**(a_DI['763m_delta']/-2.5)
    a_DI['H_contr'] = 10**(a_DI['H_delta']/-2.5)
    sz_di = 20
    alpha_di = 0.7

if include_DI_H:
    plt.scatter(a_DI['Rho(as)'],a_DI['H_contr'],color=c_h, edgecolor='k', \
        alpha=alpha_di, marker='s', s=sz_di, zorder=2, label='Known DI, H-band')

if include_DI_750_extrap:
    plt.scatter(a_DI['Rho(as)'],a_DI['763m_contr'],color=c_750, edgecolor='k', \
        marker='d', alpha=alpha_di, s=sz_di-5, zorder=2, label='Known DI, 750nm extrap.')
    if not include_DI_550_extrap:
        for ct, rho in enumerate(a_DI['Rho(as)']):
            plt.plot([rho,rho], [a_DI[ct]['763m_contr'], a_DI[ct]['H_contr']], \
            color='lightgray', linewidth=1, linestyle=':', zorder=1)

if include_DI_550_extrap:
    for ct, rho in enumerate(a_DI['Rho(as)']):
        plt.plot([rho,rho], [a_DI[ct]['547m_contr'], a_DI[ct]['H_contr']], \
            color='lightgray', linewidth=1, linestyle=':', zorder=1)
    plt.scatter(a_DI['Rho(as)'],a_DI['547m_contr'],color=c_550, edgecolor='k', \
        marker='o', alpha=alpha_di, s=sz_di-5, zorder=2, label='Known DI, 550nm extrap.')


#########################################################################
### Specific planetary systems
ProximaCenb=np.array((0.035,4E-8))
plt.plot(ProximaCenb[0],ProximaCenb[1],marker='o',color=c_h,markersize=markersize_points) ##Add proxima
plt.text(ProximaCenb[0],ProximaCenb[1],'  Proxima Cen b',color='black',horizontalalignment='left',verticalalignment='center',fontsize=10)



#########################################################################
###Add Earth and Jupiter
Earth=np.array((0.09,1*10**-10))
Jupiter=np.array((0.51, 10**-9))
plt.plot(Earth[0],Earth[1],marker='o',color='blue',markersize=markersize_points)
plt.plot(Jupiter[0],Jupiter[1],marker='o',color='#FF7800',markersize=markersize_points)
plt.text(Earth[0],Earth[1],' Earth ',color='blue',horizontalalignment='left',verticalalignment='center',fontsize=10)
plt.text(Jupiter[0],Jupiter[1],' Jupiter ',color='#FF7800',horizontalalignment='left',verticalalignment='center',fontsize=10)

#########################################################################
###Add RV planets
a_RV = np.loadtxt(path+'RVtable.txt')
arcsec_RV=a_RV[:,0]
contrast_RV=a_RV[:,1]
plt.plot(arcsec_RV,contrast_RV, 'ko', markersize=markersize_points, \
    label='Known RV, \nextrap. visible reflected')




#########################################################################
###Plot axes, tick mark ajdusting, legend, etc.
x_ticklabels=np.array((0.1,1))##arcsec
x_ticks=x_ticklabels

###To include minor tick labels include the following (note this is not an elegant solution...)
x_ticks_minor=np.array((0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0))
x_ticklabels_minor=np.array((0.03,'','','','','','','','','','','','','','','','','',''))

#y_ticks=np.array((10**-11,10**-10,10**-9,10**-8,10**-7,10**-6,10**-5,10**-4,10**-3))
#y_ticklabels=np.array(('$10^{-11}$','$10^{-10}$','$10^{-9}$','$10^{-8}$','$10^{-7}$','$10^{-6}$','$10^{-5}$','$10^{-4}$', '$10^{-3}$'))


ax1.set_ylabel('Flux ratio to host star',fontsize=14)
ax1.set_xlabel('Separation [arcsec]',fontsize=14)

ax1.set_ylim(10**-11,10**-2.9)
ax1.set_xlim(0.03,4)
ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.set_xticks(x_ticks)
ax1.set_xticklabels(x_ticklabels)
#ax1.set_yticks(y_ticks)
#ax1.set_yticklabels(y_ticklabels)

###Add grid lines
plt.grid(b=True, which='major', color='tan', linestyle='-', alpha=0.1)

####Add custom legend (probably a better way to do this?)
plt.text(0.037,2*10**-11,'Solar system planets as seen from 10 pc',color='black',horizontalalignment='left',verticalalignment='center',fontsize=8)

plt.legend(fontsize=8)

###Add minor tick labels
#ax1.set_yticks(y_ticks)
#ax1.set_yticklabels(y_ticklabels)

###Again, not an elegant solution, but it works.
ax1.set_xticks(x_ticks_minor, minor = True)
ax1.set_xticklabels(x_ticklabels_minor, minor=True)



print caption
plt.savefig(path+'flux_ratio_plot'+file_name_end+'.pdf')
#plt.savefig(path+'flux_ratio_plot'+file_name_end+'.jpg')
