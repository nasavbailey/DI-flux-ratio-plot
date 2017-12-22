#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: meshkat (original; July 31, 2017)
November-December 2017: VPB reorganized & updated detection limit curves
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rcParams
from astropy.io import ascii
import reflected_light_planets as rlp
from astropy import units as u

### User-defined options

# Which plot elements to include?
include_ELT     = False
include_HABEX   = False
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
include_BTR_spec = True # spectroscopy BTR
include_BTR_disk_to_img = True # disk mask BTR
include_special_systems = True # include, Tau Ceti and Solar System?

is_draft = True # print DRAFT on the plot?

 # colorcode detection limit curve lines by wavelength?
color_by_lambda = 'simple' # full, simple, none

###Define path where to find data and where to save plot
path = './' # leave blank or set to './' if this script is in the same folder as the data (default)
datapath = path+'data/'

########################################################################
### Plot setup

rcParams.update({'figure.autolayout': True})
#rcParams.update({'font.family':'Times New Roman'})
rcParams.update({'font.size': 12})
rcParams['mathtext.fontset'] = 'stixsans' #'stix'
rcParams['lines.solid_capstyle'] = 'butt' #don't increase line length when increasing width
rcParams['patch.linewidth'] = 0.5  # make marker edge linewidths narrower for scatter

fig=plt.figure(figsize=[6,5])#
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()

ylim = np.array([1E-11, 1E-3])
xlim = np.array([0.07,5])

rv_markersize = 35
di_markersize = 15
ccfs = 8 # instrument curve font size
ccc = 'k' # default instrument curve color
cclw = 1 # default instrument curve line width
c_pl = 'c' # color for special planetary systems (Solar System, Prox Cen, etc.)

# text about detection limit curves
ax1.text(xlim[0], ylim[0]*1.1, ' Instrument curves are \n 5$\mathdefault{\sigma}$ post-processed detection limits.',\
    color='k',horizontalalignment='left', verticalalignment='bottom',fontsize=ccfs-1)

if color_by_lambda.lower() == 'full':
    c_v = 'dodgerblue'
    c_bbvis = 'cadetblue'
    c_band3 = 'goldenrod'
    c_band4 = 'orange'
    c_yjh = 'coral'
    c_k = 'firebrick'
    c_h   = 'red'
    #c_l   = 'darkred'

    ax2.plot([1,1],[1,1],color=c_v,linewidth=cclw+2, label='< 650 nm')
    ax2.plot([1,1],[1,1],color=c_bbvis,linewidth=cclw+2, label='broadband\nvisible')
    ax2.plot([1,1],[1,1],color=c_band3,linewidth=cclw+2, label='CGI Band 3')
    ax2.plot([1,1],[1,1],color=c_band4,linewidth=cclw+2, label='CGI Band 4')
    ax2.plot([1,1],[1,1],color=c_yjh,linewidth=cclw+2, label='YJH-band')
    ax2.plot([1,1],[1,1],color=c_h,linewidth=cclw+2, label='H-band')
    ax2.plot([1,1],[1,1],color=c_k,linewidth=cclw+2, label='K-band')

elif color_by_lambda.lower() == 'simple':
    c_v = 'dodgerblue'
    c_bbvis = c_v
    c_band3 = 'orange'
    c_band4 = 'tomato'
    c_yjh = 'firebrick'
    c_h = c_yjh
    c_k = c_yjh
    #c_l = 'darkred'

    ax2.plot([1,1],[1,1],color=c_v,linewidth=cclw+2, label='< 650 nm')
    ax2.plot([1,1],[1,1],color=c_band3,linewidth=cclw+2, label='650 - 800nm')
    ax2.plot([1,1],[1,1],color=c_band4,linewidth=cclw+2, label='800 - 1000nm')
    ax2.plot([1,1],[1,1],color=c_h,linewidth=cclw+2, label='> 1000 nm')

elif color_by_lambda.lower() == 'none':
    c_v = ccc
    c_bbvis = ccc
    c_band3 = ccc
    c_band4 = ccc
    c_yjh = ccc
    c_k = ccc
    c_h   = ccc
    #c_l   = ccc

else:
    raise Exception(color_by_lambda+' is not a valid option for color_by_lambda (full/simple/none)')


########################################################################
# auto-generated caption. See README for how to comment datafiles.

caption = '** This short caption is auto-generated. DO NOT EDIT. **\n' + \
        'Please see individual datafiles for full descriptions. \n\n'

if color_by_lambda.lower() != 'none':
    caption += 'Lines and points are color coded by wavelength of observation.\n\n'

def extract_short_caption(filename):
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    for l in lines:
        if '#short caption:' in l.lower():
            return '-- '+l.split('caption:')[1].strip()+'\n\n'
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
    ax1.text(0.1,2.5*10**-7,'Future ELTs (NIR)?',color='Red',horizontalalignment='left',\
    	verticalalignment='top',rotation=-19,fontsize=ccfs)


#########################################################################
### HabEx "goal" detection limit

if include_HABEX:
    ax1.plot([0.06, 1.6],[5E-11, 5E-11],color=c_bbvis,linestyle='--',linewidth=cclw,label='')
    ax1.text(1.6,6E-11,'HabEx goal',color=c_bbvis,horizontalalignment='right',fontsize=ccfs)
    caption += '-- HabEx: Goal 5-sigma post-processed contrast.  '+\
                'IWA ~ 2.5 lambda/D @ 450nm; OWA ~ 32 l/D @ 1micron '+\
                '(source: B. Mennesson, personal communication)\n\n'


#########################################################################
### NIRCAM detection limit

if include_NIRCAM:
    fname = datapath+'jwst_nircam.txt'
    a_JWST = ascii.read(fname)
    ax1.plot(a_JWST['Rho(as)'],a_JWST['210_contr'],color=c_k,linewidth=cclw,linestyle='--',label='')
    if include_SPHERE:
        xy=[1.6, 5E-7]
        ax1.text(xy[0],xy[1], 'JWST NIRCam', color=c_k, rotation=-25, fontsize=ccfs)
        ax1.plot([1.45,xy[0]],[3E-7,xy[1]-1E-7],'k', linewidth=0.5)
    else:
        ax1.text(2,1E-7,'JWST NIRCam',color=c_k,\
            horizontalalignment='left',rotation=-30,fontsize=ccfs)

    caption += extract_short_caption(fname)

#########################################################################
### NICMOS detection limit
if include_NICMOS:
    fname = datapath+'HST_NICMOS_Min.txt' #path+'HST_NICMOS_Median.txt'
    a_NICMOS = ascii.read(fname)
    ax1.plot(a_NICMOS['Rho(as)'],a_NICMOS['F160W_contr'],color=c_h,\
        linewidth=cclw,label='')
    ax1.text(max(a_NICMOS['Rho(as)']*1.1), min(a_NICMOS['F160W_contr']), 'HST NICMOS',\
        color=c_h,horizontalalignment='right', verticalalignment='bottom', \
        rotation=-20,fontsize=ccfs)
    caption += extract_short_caption(fname)

#########################################################################
### STIS Bar5 detection limit
if include_STIS:
    fname = datapath+'HST_STIS.txt'
    a_STIS = ascii.read(fname)
    ax1.plot(a_STIS['Rho(as)'],a_STIS['KLIP_Contr'],color=c_bbvis,\
        linewidth=cclw,label='')
    ax1.text(0.2,5*10**-5,'HST STIS',color=c_bbvis,horizontalalignment='left',rotation=-40,fontsize=ccfs)
    caption += extract_short_caption(fname)

#########################################################################
### ACS detection limit

if include_ACS:
    fname = datapath+'HST_ACS.txt'
    a_ACS = ascii.read(fname)
    ax1.plot(a_ACS['Rho(as)'],a_ACS['F606W_contr'],color=c_v,linewidth=cclw,label='')
    ax1.text(4,6*10**-9,'HST ACS',color=c_v,horizontalalignment='right',rotation=-35,fontsize=ccfs)
    caption += extract_short_caption(fname)


#########################################################################
### SPHERE detection limit
if include_SPHERE:
    fname = datapath+'SPHERE_Vigan.txt'
    a_SPHERE = ascii.read(fname)
    a_SPHERE['Contrast'] = 10**(-0.4*a_SPHERE['delta'])

    # manually split into IFS and IRDIS, at 0.7", as per documentation.
    idx_yjh = a_SPHERE['Rho(as)'] <= 0.7 # IFS YJH
    idx_k12 = a_SPHERE['Rho(as)'] >= 0.7  # IRDIS K1-K2
    ax1.plot(a_SPHERE['Rho(as)'][idx_yjh], a_SPHERE['Contrast'][idx_yjh], \
        color=c_yjh, linewidth=cclw, label='')
    ax1.plot(a_SPHERE['Rho(as)'][idx_k12], a_SPHERE['Contrast'][idx_k12], \
        color=c_k, linewidth=cclw, label='')
    ax1.text(0.2,1E-6,'VLT SPHERE',color=c_k,horizontalalignment='right',fontsize=ccfs)
    ax1.text(0.14,5*10**-7,'IFS /',color=c_yjh,horizontalalignment='right',fontsize=ccfs)
    ax1.text(0.14,5*10**-7,' IRDIS',color=c_k,horizontalalignment='left',fontsize=ccfs)
    caption += extract_short_caption(fname)

#########################################################################
### GPI H-band
if include_GPI:
    fname = datapath+'GPIES_T-type_contrast_curve_2per.txt'
    a_GPI = ascii.read(fname)
    ax1.plot(a_GPI['Rho(as)'],a_GPI['H_contr_60min'],color=c_h,linewidth=cclw,label='')
    ax1.text(0.17,0.8*10**-5,'Gemini GPI',color=c_h,horizontalalignment='left',rotation=-24,fontsize=ccfs)
    caption += extract_short_caption(fname)

#########################################################################
### WFIRST
a_e = np.loadtxt(datapath+'exoplanet_mode.csv',delimiter=',')
a_d = np.loadtxt(datapath+'disk_mode.csv',delimiter=',')


arcsec_exoplanet=a_e[:,1]
contrast_exoplanet=a_e[:,2]
arcsec_disk=a_d[:,1]
contrast_disk=a_d[:,2]

ax1.plot(arcsec_exoplanet,contrast_exoplanet,color='m', linewidth=cclw+4)
ax1.plot(arcsec_disk,contrast_disk, color='m', linewidth=cclw+4)
caption += '-- WFIRST curves are pre-WEITR L3 requirements for 5-sigma, post-processed detection limits.\n\n'
ax1.text(1.3, 2E-9, 'WFIRST\nCGI\npre-WEITR', color='m', horizontalalignment='left',\
    verticalalignment='center', fontsize=ccfs+1, weight='bold')

######Add Technical requirement line and text
if include_BTR_img:
    fname = datapath+'WFIRST_BTR_imaging.txt'
    btr_img = ascii.read(fname)
    ax1.plot(btr_img['Rho(as)'], btr_img['Band1_contr_snr5'], color=c_v, linewidth=cclw+4, label='')
    ax1.text(btr_img['Rho(as)'][0], btr_img['Band1_contr_snr5'][0], 'img \n BTR ', color=c_v,\
        horizontalalignment='right', verticalalignment='center', weight='bold', fontsize=ccfs+1)
    caption += extract_short_caption(fname)

if include_BTR_disk_to_img:
    fname = datapath+'WFIRST_BTR_disk.txt'
    btr_disk = ascii.read(fname)
    ax1.plot(btr_disk['Rho(as)'], btr_disk['Band4_pt_contr_snr5'], color=c_band4, linewidth=cclw+4, label='')
    ax1.text(0.9*btr_disk['Rho(as)'][-1], 1.1*btr_disk['Band4_pt_contr_snr5'][-1], 'disk\n BTR', color=c_band4,\
        horizontalalignment='right', verticalalignment='center', weight='bold', fontsize=ccfs+1)
    caption += '-- WFIRST extended source (disk) BTR translated to point source flux ratio\n'
    caption += extract_short_caption(fname)

if include_BTR_spec:
    import matplotlib.patheffects as path_effects
    fname = datapath+'WFIRST_BTR_spec.txt'
    btr_img = ascii.read(fname)
    if include_BTR_disk_to_img:  # draw a shadow under the line to make it easier to see if overlapping other BTR
        p = [path_effects.SimpleLineShadow(offset=(1, -1)), path_effects.Normal()]
    else:
        p = [path_effects.Normal()]
    ax1.plot(btr_img['Rho(as)'], btr_img['Band3_contr_snr5'], color=c_band3, \
        linewidth=cclw+2, label='', path_effects=p)
    ax1.text(btr_img['Rho(as)'][0], btr_img['Band3_contr_snr5'][0], 'spec \n BTR ', color=c_band3,\
        horizontalalignment='left', verticalalignment='bottom', weight='bold', fontsize=ccfs+1)
    caption += extract_short_caption(fname)


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
    alpha_di = 0.7
    caption += extract_short_caption(fname)

if include_DI_H:
    ax1.scatter(a_DI['Rho(as)'],a_DI['H_contr'],color=c_h, edgecolor='k', \
        alpha=alpha_di, marker='s', s=di_markersize, zorder=2, label='DI, 1.6 $\mathdefault{\mu} $m')

if include_DI_750_extrap:
    ax1.scatter(a_DI['Rho(as)'],a_DI['763m_contr'],color=c_band3, edgecolor='k', \
        marker='d', alpha=alpha_di, s=di_markersize+15, zorder=2, label='DI, 750nm pred.')
    if not include_DI_550_extrap:
        for ct, rho in enumerate(a_DI['Rho(as)']):
            ax1.plot([rho,rho], [a_DI[ct]['763m_contr'], a_DI[ct]['H_contr']], \
            color='lightgray', linewidth=1, linestyle=':', zorder=1)

if include_DI_550_extrap:
    for ct, rho in enumerate(a_DI['Rho(as)']):
        ax1.plot([rho,rho], [a_DI[ct]['547m_contr'], a_DI[ct]['H_contr']], \
            color='lightgray', linewidth=1, linestyle=':', zorder=1)
    ax1.scatter(a_DI['Rho(as)'],a_DI['547m_contr'],color=c_v, edgecolor='k', \
        marker='o', alpha=alpha_di, s=di_markersize+15, zorder=2, label='DI, 550nm pred.')


#########################################################################
### Specific planetary systems

if include_special_systems:
    # Proxima Centari b
    #sma = 0.05*u.au
    #flux_ratio = rlp.calc_lambert_flux_ratio(sma=sma, rp=1.3**(1./3)*u.earthRad,\
    #    orb_ang=0*u.degree,albedo=0.3, inclin=0*u.degree)
    #rho = (sma/1.3*u.pc).value
    #ax1.plot(rho,flux_ratio,marker='^', color=c_pl)
    #ax1.text(rho,flux_ratio,'  Proxima Cen b',color=c_pl,\
    #    horizontalalignment='left',verticalalignment='center',fontsize=ccfs)
    #caption += '-- Prox Cen b: At quadrature, albedo = 0.3, radius = (M/Me)^(1/3) * Re, circular orbit.\n\n'


    # Tau Ceti
    tc_dist = 3.65*u.pc
    albedo = 0.35
    # make a separate upper x axis for physical separation of this system
    ax3 = ax1.twiny()
    ax3.set_ylim(ylim)
    ax3.set_xlim(xlim * tc_dist.value)
    ax3.set_xscale('log')
    ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    ax3.set_xlabel('Semi-major axis for Tau Ceti (au)')

    # Tau Ceti f
    sma = 1.334*u.au
    flux_ratio = rlp.calc_lambert_flux_ratio(sma=sma, rp=3.9**(1./3)*u.earthRad,\
        orb_ang=0*u.degree,albedo=albedo, inclin=0*u.degree)
    rho = (sma/tc_dist).value
    ax3.scatter(sma, flux_ratio,marker='^', color=c_pl, edgecolor='k')
    ax3.text(sma.value,flux_ratio,'  Tau Ceti f',color=c_pl,\
        horizontalalignment='left',verticalalignment='center',fontsize=ccfs)

    # Tau Ceti e
    sma = 0.538*u.au
    flux_ratio = rlp.calc_lambert_flux_ratio(sma=sma, rp=3.9**(1./3)*u.earthRad,\
        orb_ang=0*u.degree,albedo=albedo, inclin=0*u.degree)
    rho = (sma/tc_dist).value
    #ax1.plot(rho,flux_ratio,marker='^', color=c_pl)
    ax3.scatter(sma, flux_ratio,marker='^', color=c_pl, edgecolor='k')
    ax3.text(sma.value,flux_ratio,'Tau Ceti e  ',color=c_pl,\
        horizontalalignment='right',verticalalignment='center',fontsize=ccfs)

    caption += '-- Tau Ceti e&f. At quadrature, albedo = ' + str(albedo) +\
        ', radius = (M/Me)^(1/3) * Re, circular orbits.\n\n'

    # Earth & Jupiter
    #earthRatio = rlp.calc_lambert_flux_ratio(sma=1.*u.au, rp=1.*u.earthRad,\
    #    orb_ang=0*u.degree,albedo=0.367, inclin=0*u.degree)
    earthRatio = 1.0E-10 # use the "standard" value
    jupiterRatio = rlp.calc_lambert_flux_ratio(sma=5.*u.au, rp=1.*u.jupiterRad,\
        orb_ang=0*u.degree,albedo=0.52, inclin=0*u.degree)
    caption += '-- Earth (==1E-10) & Jupiter (albedo=0.52) at quadrature as seen from 10 pc. '+\
                '(Jupiter albedo: Traub & Oppenheimer, '+\
                'Direct Imaging chapter of Seager Exoplanets textbook, Table 3)\n\n'

    ax1.scatter(0.1,earthRatio,marker='$\\bigoplus$',color=c_pl, s=rv_markersize)
    ax1.text(0.1,earthRatio,'  Earth ',color=c_pl,\
        horizontalalignment='left',verticalalignment='center',fontsize=ccfs)

    ax1.scatter(0.5,jupiterRatio,marker='v',color=c_pl,  edgecolor='k', s=rv_markersize)
    ax1.text(0.5,jupiterRatio,'  Jupiter ',color=c_pl,\
    horizontalalignment='left',verticalalignment='center',fontsize=ccfs)

    ax1.text(xlim[1], ylim[0]*1.1, 'Solar System as seen from 10pc. ',\
    color='k',horizontalalignment='right', verticalalignment='bottom',fontsize=ccfs-1)


#########################################################################
###Add RV planets
fname = datapath+'reflected_light_table.txt'
tmp = ascii.read(fname)
idx_rv = tmp['pl_discmethod'] == "Radial Velocity"
ax1.scatter(tmp[idx_rv]['sma_arcsec'],tmp[idx_rv]['Fp/F*_quad'],s=rv_markersize,\
    color='dimgray', edgecolor='k', marker='^', label='RV, reflected light', zorder=5)
caption += extract_short_caption(fname)


#########################################################################
###Plot axes, tick mark ajdusting, legend, etc.

if is_draft:
    from datetime import date
    ax1.text(xlim[1], ylim[0]*2, "DRAFT  "+str(date.today()) + ' ', \
        horizontalalignment='right',verticalalignment='bottom',color='red',weight='bold')

ax1.legend(fontsize=7, loc='upper right')

ax1.grid(b=True, which='major', color='tan', linestyle='-', alpha=0.1)

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xticks([0.05, 0.1, 0.5, 1, 5])
ax1.set_ylim(ylim)
ax1.set_xlim(xlim)

# write x axis in scalar notation instead of powers of 10
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1g'))

if color_by_lambda.lower() != 'none':
    second_legend = ax2.legend(loc='upper left', fontsize=7, title='Bandpass')
    second_legend.get_title().set_fontsize(8)
    ax2.set_ylim(ylim)
    ax2.set_xlim(xlim)
    #ax2.set_yscale('log')
    #ax2.set_xscale('log')
    ax2.set_yticklabels([])
    ax2.yaxis.set_ticks_position('none')

ax1.set_ylabel('Flux ratio to host star')
ax1.set_xlabel('Separation [arcsec]')



plt.tight_layout()

f = open(path+'auto_caption.txt','w') #open and overwrite existing
f.write(caption)
f.close()
plt.savefig(path+'flux_ratio_plot.pdf')
#plt.savefig(path+'flux_ratio_plot.jpg')
