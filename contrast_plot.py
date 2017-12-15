#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: meshkat (original; July 31, 2017)
November-December 2017: VPB reorganized & updated contrast curves
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rcParams
from astropy.io import ascii
import reflected_light_planets as rlp
from astropy import units as u

### User-defined options

# Which contrast curves to include?
include_ELT     = False
include_HABEX   = True
include_ACS     = True
include_NICMOS  = True
include_STIS    = True
include_NIRCAM  = True
include_SPHERE  = True
include_GPI     = True
include_DI_H    = True # real H-band contrasts of known directly imaged planets
include_DI_750_extrap = True # COND/BT-Settl model extrapolations to ~750nm
include_DI_550_extrap = True # COND/BT-Settl model extrapolations to ~550nm
include_BTR_img = True # imaging BTR
include_BTR_disk_to_img = True # disk mask BTR

is_draft = True # print DRAFT on the plot?

 # colorcode contrast curve lines by wavelength?
color_by_lambda = 'simple' # full, simple, none

###Define path where to find data and where to save plot
path = './' # leave blank or set to './' if this script is in the same folder as the data (default)
datapath = path+'data/'

########################################################################
### Plot setup

rcParams.update({'figure.autolayout': True})
#rcParams.update({'font.family':'Times New Roman'})
rcParams.update({'font.size': 12})
rcParams['mathtext.fontset'] = 'stix'
rcParams['lines.solid_capstyle'] = 'butt' #don't increase line length when increasing width

fig=plt.figure()#figsize=[7,4.5]
ax1=fig.add_subplot(111)
#ax2 = ax1.twinxy()

markersize_points=4
ccfs = 8 # contrast curve font size
ccc = 'darkviolet' # default contrast curve color
cclw = 1.5 # default contrast curve line width
c_pl = 'c' # color for special planetary systems (Solar System, Prox Cen, etc.)

if color_by_lambda.lower() == 'full':
    c_v = 'dodgerblue'
    c_bbvis = 'cadetblue'
    c_750 = 'goldenrod'
    c_yjh = 'coral'
    c_k = 'firebrick'
    c_h   = 'red'
    #c_l   = 'darkred'

    line1, = plt.plot([1,1],[1,1],color=c_v,linewidth=cclw, label='V/550nm')
    line2, = plt.plot([1,1],[1,1],color=c_bbvis,linewidth=cclw, label='broadband visible')
    line3, = plt.plot([1,1],[1,1],color=c_750,linewidth=cclw, label='750nm')
    line4, = plt.plot([1,1],[1,1],color=c_yjh,linewidth=cclw, label='YJH-band')
    line5, = plt.plot([1,1],[1,1],color=c_h,linewidth=cclw, label='H-band')
    line6, = plt.plot([1,1],[1,1],color=c_k,linewidth=cclw, label='K-band')

    first_legend = plt.legend(handles=[line1, line2, line3, line4, line5, line6], \
        loc='upper right', fontsize=7, title='Bandpass')
    first_legend.get_title().set_fontsize('8')

elif color_by_lambda.lower() == 'simple':
    c_v = 'dodgerblue'
    c_bbvis = c_v
    c_750 = 'goldenrod'
    c_yjh = 'firebrick'
    c_h = c_yjh
    c_k = c_yjh
    #c_l = 'darkred'

    line1, = plt.plot([1,1],[1,1],color=c_v,linewidth=cclw, label='Blue optical')
    line2, = plt.plot([1,1],[1,1],color=c_750,linewidth=cclw, label='Red optical')
    line3, = plt.plot([1,1],[1,1],color=c_h,linewidth=cclw, label='near IR')

    first_legend = plt.legend(handles=[line1, line2, line3], \
        loc='upper right', fontsize=7, title='Bandpass')
    first_legend.get_title().set_fontsize('8')


elif color_by_lambda.lower() == 'none':
    c_v = ccc
    c_750 = ccc
    c_yjh = ccc
    c_k = ccc
    c_h   = ccc
    #c_l   = ccc

else:
    raise Exception(color_by_lambda+' is not a valid option for color_by_lambda (full/simple/none)')

#ax = plt.gca().add_artist(first_legend)


########################################################################
# auto-generated caption. See README for how to comment datafiles.

caption = '** This short caption is auto-generated. DO NOT EDIT. **\n' + \
        'Please see individual datafiles for full descriptions. \n\n'

def extract_short_caption(filename):
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    for l in lines:
        if '#short caption:' in l.lower():
            return '-- '+l.split('caption:')[1].strip()+'\n'
    # if no caption in text file
    print '\n**** WARNING **** no caption for '+filename+'\n'
    return ''



########################################################################
### ELT

if include_ELT:
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
    plt.plot([0.06, 1.6],[5E-11, 5E-11],color=c_bbvis,linestyle='--',linewidth=cclw-0.5,label='')
    plt.text(1.6,6E-11,'HabEx',color=c_bbvis,horizontalalignment='right',fontsize=ccfs)
    caption += '-- HabEx: Goal 5-sigma post-processed contrast.  '+\
                'IWA ~ 2.5 lambda/D @ 450nm; OWA ~ 32 l/D @ 1micron '+\
                '(source: B. Mennesson, personal communication)\n'


#########################################################################
### NIRCAM contrast curve

if include_NIRCAM:
    fname = datapath+'jwst_nircam.txt'
    a_JWST = ascii.read(fname)
    plt.plot(a_JWST['Rho(as)'],a_JWST['210_contr'],color=c_k,linewidth=cclw-0.5,linestyle='--',label='')
    if include_SPHERE:
        xy=[1.6, 5E-7]
        plt.text(xy[0],xy[1], 'JWST NIRCam', color=c_k, rotation=-25, fontsize=ccfs)
        plt.plot([1.45,xy[0]],[3E-7,xy[1]-1E-7],'k', linewidth=0.5)
    else:
        plt.text(2,1E-7,'JWST NIRCam',color=c_k,\
            horizontalalignment='left',rotation=-30,fontsize=ccfs)

    caption += extract_short_caption(fname)

#########################################################################
### NICMOS contrast curve
if include_NICMOS:
    fname = datapath+'HST_NICMOS_Min.txt' #path+'HST_NICMOS_Median.txt'
    a_NICMOS = ascii.read(fname)
    plt.plot(a_NICMOS['Rho(as)'],a_NICMOS['F160W_contr'],color=c_h,\
        linewidth=cclw,label='')
    plt.text(max(a_NICMOS['Rho(as)']*1.1),3.5E-6,'HST NICMOS',\
        color=c_h,horizontalalignment='right',rotation=-20,fontsize=ccfs)
    caption += extract_short_caption(fname)

#########################################################################
### STIS Bar5 contrast curve
if include_STIS:
    fname = datapath+'HST_STIS.txt'
    a_STIS = ascii.read(fname)
    plt.plot(a_STIS['Rho(as)'],a_STIS['KLIP_Contr'],color=c_bbvis,\
        linewidth=cclw,label='')
    plt.text(0.2,5*10**-5,'HST STIS',color=c_bbvis,horizontalalignment='left',rotation=-40,fontsize=ccfs)
    caption += extract_short_caption(fname)

#########################################################################
### ACS contrast curve

if include_ACS:
    fname = datapath+'HST_ACS.txt'
    a_ACS = ascii.read(fname)
    plt.plot(a_ACS['Rho(as)'],a_ACS['F606W_contr'],color=c_v,linewidth=cclw,label='')
    plt.text(3.8,6*10**-9,'HST ACS',color=c_v,horizontalalignment='right',rotation=-35,fontsize=ccfs)
    caption += extract_short_caption(fname)


#########################################################################
### SPHERE contrast curve
if include_SPHERE:
    fname = datapath+'SPHERE_Vigan.txt'
    a_SPHERE = ascii.read(fname)
    a_SPHERE['Contrast'] = 10**(-0.4*a_SPHERE['delta'])

    # manually split into IFS and IRDIS, at 0.7", as per documentation.
    idx_yjh = a_SPHERE['Rho(as)'] <= 0.7 # IFS YJH
    idx_k12 = a_SPHERE['Rho(as)'] >= 0.7  # IRDIS K1-K2
    plt.plot(a_SPHERE['Rho(as)'][idx_yjh], a_SPHERE['Contrast'][idx_yjh], \
        color=c_yjh, linewidth=cclw, label='')
    plt.plot(a_SPHERE['Rho(as)'][idx_k12], a_SPHERE['Contrast'][idx_k12], \
        color=c_k, linewidth=cclw, label='')
    plt.text(0.2,1E-6,'VLT SPHERE',color=c_k,horizontalalignment='right',fontsize=ccfs)
    plt.text(0.14,5*10**-7,'IFS /',color=c_yjh,horizontalalignment='right',fontsize=ccfs)
    plt.text(0.14,5*10**-7,' IRDIS',color=c_k,horizontalalignment='left',fontsize=ccfs)
    caption += extract_short_caption(fname)

#########################################################################
### GPI H-band
if include_GPI:
    fname = datapath+'GPIES_T-type_contrast_curve_2per.txt'
    a_GPI = ascii.read(fname)
    plt.plot(a_GPI['Rho(as)'],a_GPI['H_contr'],color=c_h,linewidth=cclw,label='')
    plt.text(0.17,1*10**-5,'Gemini GPI',color=c_h,horizontalalignment='left',rotation=-20,fontsize=ccfs)
    caption += extract_short_caption(fname)

#########################################################################
### WFIRST
a_e = np.loadtxt(datapath+'exoplanet_mode.csv',delimiter=',')
a_d = np.loadtxt(datapath+'disk_mode.csv',delimiter=',')


arcsec_exoplanet=a_e[:,1]
contrast_exoplanet=a_e[:,2]
arcsec_disk=a_d[:,1]
contrast_disk=a_d[:,2]

plt.plot(arcsec_exoplanet,contrast_exoplanet,color='black',linestyle='--',\
    linewidth=cclw+0.5,label='WFIRST')
plt.plot(arcsec_disk,contrast_disk,color='black',linestyle='--',\
    linewidth=cclw+0.5)
caption += '-- WFIRST lines are pre-WEITR L3 requirements for 5-sigma, post-processed contrast.\n'

######Add Technical requirement line and text
if include_BTR_img:
    plt.plot([0.23, 0.4], [0.5*5E-8, 0.5*5E-8], color=c_v, linewidth=cclw+2, alpha=0.7)
    plt.text(0.23,2.5*10**-8,'BTR1 ',color=c_v,horizontalalignment='right',\
        weight='bold',fontsize=ccfs)
    caption += '-- BTR1: imaging BTR.\n'

if include_BTR_disk_to_img:
    plt.plot([0.25, 0.95], [0.5*5E-8, 0.5*5E-8], color=c_750, linewidth=cclw+4, alpha=0.6)
    plt.text(0.95,2.5*10**-8,' BTR3',color=c_750,horizontalalignment='left',\
        weight='bold',fontsize=ccfs)
    caption += '-- BTR3: extended object sensitivity BTR translate to point source sensitivity.\n'



#########################################################################
###### -------------------- planets -------------------------  ##########
#########################################################################

#########################################################################
### Self luminous directly imaged planets

if include_DI_H or include_DI_extrap:
    fname = datapath+'DIplanets.txt'
    a_DI = ascii.read(fname)
    a_DI['547m_contr'] = 10**(a_DI['547m_delta']/-2.5)
    a_DI['763m_contr'] = 10**(a_DI['763m_delta']/-2.5)
    a_DI['H_contr'] = 10**(a_DI['H_delta']/-2.5)
    sz_di = 15
    alpha_di = 1
    caption += extract_short_caption(fname)

if include_DI_H:
    plt.scatter(a_DI['Rho(as)'],a_DI['H_contr'],color=c_h, edgecolor='k', \
        alpha=alpha_di, marker='s', s=sz_di, zorder=2, label='DI, H-band')

if include_DI_750_extrap:
    plt.scatter(a_DI['Rho(as)'],a_DI['763m_contr'],color=c_750, edgecolor='k', \
        marker='d', alpha=alpha_di, s=sz_di+15, zorder=2, label='DI, 750nm extrap.')
    if not include_DI_550_extrap:
        for ct, rho in enumerate(a_DI['Rho(as)']):
            plt.plot([rho,rho], [a_DI[ct]['763m_contr'], a_DI[ct]['H_contr']], \
            color='lightgray', linewidth=1, linestyle=':', zorder=1)

if include_DI_550_extrap:
    for ct, rho in enumerate(a_DI['Rho(as)']):
        plt.plot([rho,rho], [a_DI[ct]['547m_contr'], a_DI[ct]['H_contr']], \
            color='lightgray', linewidth=1, linestyle=':', zorder=1)
    plt.scatter(a_DI['Rho(as)'],a_DI['547m_contr'],color=c_v, edgecolor='k', \
        marker='o', alpha=alpha_di, s=sz_di+15, zorder=2, label='DI, 550nm extrap.')


#########################################################################
### Specific planetary systems
# Proxima Centari b
sma = 0.05*u.au
flux_ratio = rlp.calc_lambert_flux_ratio(sma=sma, rp=1.3**(1./3)*u.earthRad,\
    orb_ang=0*u.degree,albedo=0.3, inclin=0*u.degree)
rho = (sma/1.3*u.pc).value
plt.plot(rho,flux_ratio,marker='+', color=c_pl)
plt.text(rho,flux_ratio,'  Proxima Cen b',color=c_pl,\
    horizontalalignment='left',verticalalignment='center',fontsize=ccfs)

# Tau Ceti f
sma = 1.334*u.au
flux_ratio = rlp.calc_lambert_flux_ratio(sma=sma, rp=3.9**(1./3)*u.earthRad,\
    orb_ang=0*u.degree,albedo=0.3, inclin=0*u.degree)
rho = (sma/3.65*u.pc).value
plt.plot(rho,flux_ratio,marker='+', color=c_pl)
plt.text(rho,flux_ratio,'  Tau Ceti f',color=c_pl,\
    horizontalalignment='left',verticalalignment='center',fontsize=ccfs)

# Tau Ceti e
sma = 0.538*u.au
flux_ratio = rlp.calc_lambert_flux_ratio(sma=sma, rp=3.9**(1./3)*u.earthRad,\
    orb_ang=0*u.degree,albedo=0.3, inclin=0*u.degree)
rho = (sma/3.65*u.pc).value
plt.plot(rho,flux_ratio,marker='+', color=c_pl)
plt.text(rho,flux_ratio,'Tau Ceti e  ',color=c_pl,\
    horizontalalignment='right',verticalalignment='center',fontsize=ccfs)

caption += '-- Other noteable planetary systems: Prox Cen b, Tau Ceti e&f. '+\
            'At quadrature, albedo = 0.3, radius = (M/Me)^(1/3) * Re\n'

#########################################################################
###Add Earth and Jupiter
earthRatio = rlp.calc_lambert_flux_ratio(sma=1.*u.au, rp=1.*u.earthRad,\
    orb_ang=0*u.degree,albedo=0.3, inclin=0*u.degree)
jupiterRatio = rlp.calc_lambert_flux_ratio(sma=5.*u.au, rp=1.*u.jupiterRad,\
    orb_ang=0*u.degree,albedo=0.5, inclin=0*u.degree)


plt.plot(0.1,earthRatio,marker='$\\bigoplus$',color=c_pl,markersize=10)
plt.text(0.1,earthRatio,'  Earth ',color=c_pl,\
    horizontalalignment='left',verticalalignment='center',fontsize=ccfs)

plt.plot(0.5,jupiterRatio,marker='+',color=c_pl,markersize=markersize_points)
plt.text(0.5,jupiterRatio,'  Jupiter ',color=c_pl,\
horizontalalignment='left',verticalalignment='center',fontsize=ccfs)
caption += '-- Earth & Jupiter at quadrature as seen from 10 pc. '+\
            'Albedos of 0.3 and 0.5, respectively.\n'


#########################################################################
###Add RV planets
tmp = ascii.read(datapath+'reflected_light_table.txt')
idx_rv = tmp['pl_discmethod'] == "Radial Velocity"
plt.scatter(tmp[idx_rv]['sma_arcsec'],tmp[idx_rv]['Fp/F*_quad'],color='k',s=30,\
    marker='^', label='RV, extrap.\nreflected light')



#########################################################################
###Plot axes, tick mark ajdusting, legend, etc.

if is_draft:
    from datetime import date
    plt.title("DRAFT  "+str(date.today()), color='red',weight='bold')

plt.legend(fontsize=7)

plt.grid(b=True, which='major', color='tan', linestyle='-', alpha=0.1)

ax1.set_ylim(1E-11, 2E-3)
ax1.set_xlim(0.03,4)
ax1.set_yscale('log')
ax1.set_xscale('log')


ax1.set_ylabel('Flux ratio to host star')
ax1.set_xlabel('Separation [arcsec]')


###To include minor tick labels include the following (note this is not an elegant solution...)
#x_ticks_minor=np.array((0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0))
#x_ticklabels_minor=np.array((0.03,'','','','','','','','','','','','','','','','','',''))
#ax1.set_xticks(x_ticks_minor, minor = True)
#ax1.set_xticklabels(x_ticklabels_minor, minor=True)

# write x axis in scalar notation instead of powers of 10
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

plt.tight_layout()

f = open(path+'caption.txt','w') #open and overwrite existing
f.write(caption)
f.close()
plt.savefig(path+'flux_ratio_plot.pdf')
#plt.savefig(path+'flux_ratio_plot.jpg')
