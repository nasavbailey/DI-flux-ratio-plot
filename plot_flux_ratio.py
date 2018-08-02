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
import yaml
from glob import glob
from datetime import date
import os

### import YAML file of user-defined options
cfglist = glob('*yml')
if len(cfglist) > 1:
    raise Exception(('Mulitple YAML config files present: %s.\nMove all but desired config file to another folder.')%flist)
if len(cfglist) == 0:
    raise Exception('No .yml config files present. Place your .yml file in the same directory as plot_flux_ratio.py.')
with open(cfglist[0],'r') as f:
    cfg = yaml.load(f)
    cfgname = cfglist[0].split('.yml')[0]

###Define path where to find data and where to save plot
datapath = './data/'

### hardcoded constants
d_tel = 2.4 * u.m

########################################################################
### Plot setup

rcParams['figure.autolayout'] = True
rcParams['font.size'] = cfg['plot_font_size']
rcParams['mathtext.fontset'] = cfg['math_font']
rcParams['lines.solid_capstyle'] = 'butt' #don't increase line length when increasing width
rcParams['patch.linewidth'] = cfg['marker_edge_width']  # make marker edge linewidths narrower for scatter

fig=plt.figure(figsize=[cfg['fig_width'], cfg['fig_height']])#
ax1 = fig.add_subplot(111)

xlim = np.array([float(cfg['x0']), float(cfg['x1'])])
ylim = np.array([float(cfg['y0']), float(cfg['y1'])])

ccfs = cfg['label_font_size']
lw1 = cfg['other_linewidth']
lw2 = cfg['wfirst_linewidth']


if cfg['color_by_lambda'].lower() == 'full':
    c_v = 'dodgerblue'
    c_bbvis = 'cadetblue'
    c_band3 = 'goldenrod'
    c_band4 = 'orange'
    c_yjh = 'coral'
    c_k = 'firebrick'
    c_h   = 'red'
    c_pl = 'c'

    ax2 = ax1.twinx()

    if cfg['ACS'] or cfg['req_img'] or cfg['DI_B1_pred']:
        ax2.plot([1,1],[1,1],color=c_v,linewidth=lw1+2, label='< 650 nm')
    if cfg['STIS']:
        ax2.plot([1,1],[1,1],color=c_bbvis,linewidth=lw1+2, label='broadband\nvisible')
    if cfg['req_spec'] or cfg['pred_spec'] or cfg['DI_B3_pred']:
        ax2.plot([1,1],[1,1],color=c_band3,linewidth=lw1+2, label='CGI Band 3')
    if cfg['req_disk']:
        ax2.plot([1,1],[1,1],color=c_band4,linewidth=lw1+2, label='CGI Band 4')
    if cfg['SPHERE']:
        ax2.plot([1,1],[1,1],color=c_yjh,linewidth=lw1+2, label='YJH-band')
    if cfg['GPI'] or cfg['NICMOS'] or cfg['DI_H']:
        ax2.plot([1,1],[1,1],color=c_h,linewidth=lw1+2, label='H-band')
    if cfg['SPHERE'] or cfg['NIRCAM']:
        ax2.plot([1,1],[1,1],color=c_k,linewidth=lw1+2, label='K-band')

elif cfg['color_by_lambda'].lower() == 'simple':
    c_v = 'dodgerblue'
    c_bbvis = c_v
    c_band3 = 'orange'
    c_band4 = 'tomato'
    c_yjh = 'firebrick'
    c_h = c_yjh
    c_k = c_yjh
    c_pl = 'c'

    ax2 = ax1.twinx()
    if cfg['HABEX'] or cfg['ACS'] or cfg['STIS'] or cfg['DI_B1_pred']:
        ax2.plot([1,1],[1,1],color=c_v,linewidth=lw1+2, label='< 650 nm')
    if cfg['DI_B3_pred'] or cfg['req_spec'] or cfg['pred_spec']:
        ax2.plot([1,1],[1,1],color=c_band3,linewidth=lw1+2, label='650 - 800nm')
    if cfg['req_disk'] or cfg['pred_disk']:
        ax2.plot([1,1],[1,1],color=c_band4,linewidth=lw1+2, label='800 - 1000nm')
    if cfg['GPI'] or cfg['SPHERE'] or cfg['NIRCAM'] or cfg['NICMOS'] or cfg['DI_H']:
        ax2.plot([1,1],[1,1],color=c_h,linewidth=lw1+2, label='> 1000 nm')


elif cfg['color_by_lambda'].lower() == 'minimal':
    c_v = 'dodgerblue'
    c_bbvis = c_v
    c_band3 = c_v
    c_band4 = c_v
    c_yjh = 'firebrick'
    c_h = c_yjh
    c_k = c_yjh
    c_pl = 'c'

    ax2 = ax1.twinx()
    if cfg['HABEX'] or cfg['ACS'] or cfg['STIS'] or cfg['DI_B1_pred'] or \
    cfg['DI_B3_pred'] or cfg['req_spec'] or cfg['pred_spec'] or  cfg['req_disk']:
        ax2.plot([1,1],[1,1],color=c_band4,linewidth=lw1+2, label='< 1000 nm')
    if cfg['GPI'] or cfg['SPHERE'] or cfg['NIRCAM'] or cfg['NICMOS']:
        ax2.plot([1,1],[1,1],color=c_h,linewidth=lw1+2, label='> 1000 nm')


elif cfg['color_by_lambda'].lower() == 'none':
    ccc = 'k'
    c_v = ccc
    c_bbvis = ccc
    c_band3 = ccc
    c_band4 = ccc
    c_yjh = ccc
    c_k = ccc
    c_h   = ccc
    c_pl = ccc

else:
    raise Exception(cfg['color_by_lambda']+' is not a valid option for color_by_lambda (full/simple/none)')


# text about detection limit curves
if cfg['ELT'] or cfg['HABEX'] or cfg['NIRCAM'] or cfg['NICMOS'] or cfg['STIS'] \
    or cfg['ACS'] or cfg['SPHERE'] or cfg['GPI']:
    ax1.text(xlim[-1], ylim[0]*1.1, \
    ' Instrument curves are 5$\mathdefault{\sigma}$ post-processed detection limits. ',\
    horizontalalignment='right', verticalalignment='bottom',\
    fontsize=ccfs+1, color='k', weight='bold')#, backgroundcolor='white')

########################################################################
# auto-generated caption. See README for how to comment datafiles.
# auto-generated caption
caption = '** This short caption is auto-generated. DO NOT EDIT. **\n' + \
        'Please see individual datafiles for full descriptions. \n'
caption += 'This file was generated on %s\n'%str(date.today())
caption += 'Config file used = %s\n\n'%cfglist[0]

if cfg['color_by_lambda'].lower() != 'none':
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


#########################################################################
###### --------- instrument detection limits-----------------  ##########
#########################################################################

########################################################################
### ELT guess

if cfg['ELT']:
    range_x = np.array((0.03, 1))
    pessimistic_y = np.array((1E-5, 1E-8))
    optimistic_y=np.array((1E-8, 1E-9))
    ax1.plot(range_x, pessimistic_y, color='pink', linestyle='--', linewidth=lw1)
    ax1.plot(range_x, optimistic_y, color='pink', linestyle='--', linewidth=lw1)
    ax1.fill_between(range_x, pessimistic_y, optimistic_y, color='pink', alpha='0.2')

    ax1.text(0.08, 1E-8, 'ELT goal', color='coral', horizontalalignment='left',\
        verticalalignment='top', rotation=-8, fontsize=ccfs)

    caption += '-- ELT goal: Possible range of near-IR post-processed detection limits for ' + \
                'next generation extremely large telescopes. \n\n'

#########################################################################
### HabEx "goal" detection limit

if cfg['HABEX']:
    ax1.plot([0.06, 1.65],[5E-11, 5E-11],color=c_bbvis,linestyle='--',linewidth=lw1,label='')
    ax1.text(1.6,6E-11,'HabEx goal',color=c_bbvis,horizontalalignment='right',fontsize=ccfs)
    caption += '-- HabEx: Goal 5-sigma post-processed contrast.  '+\
                'IWA ~ 2.5 lambda/D @ 450nm; OWA ~ 32 l/D @ 1micron '+\
                '(source: B. Mennesson, personal communication)\n\n'


#########################################################################
### NIRCAM detection limit

if cfg['NIRCAM']:
    fname = datapath+'jwst_nircam.txt'
    a_JWST = ascii.read(fname)
    ax1.plot(a_JWST['Rho(as)'],a_JWST['210_contr'],color=c_k,linewidth=lw1,linestyle='--',label='')
    if cfg['SPHERE']:
        xy=[4, 3E-8]
        ax1.text(xy[0],xy[1], 'JWST NIRCam', color=c_k, rotation=-30, fontsize=ccfs, \
            verticalalignment='bottom', horizontalalignment='right')
        ax1.plot([xy[0]*0.95,xy[0]],[xy[1],1.1E-8],'k', linewidth=0.5)
    else:
        ax1.text(2,1E-7,'JWST NIRCam',color=c_k,\
            horizontalalignment='left',rotation=-30,fontsize=ccfs)

    caption += extract_short_caption(fname)


#########################################################################
### NICMOS detection limit

if cfg['NICMOS']:
    fname = datapath+'HST_NICMOS_Min.txt' #path+'HST_NICMOS_Median.txt'
    a_NICMOS = ascii.read(fname)
    ax1.plot(a_NICMOS['Rho(as)'],a_NICMOS['F160W_contr'],color=c_h,\
        linewidth=lw1,label='')
    ax1.text(max(a_NICMOS['Rho(as)']*1.1), min(a_NICMOS['F160W_contr']), 'HST NICMOS',\
        color=c_h,horizontalalignment='right', verticalalignment='bottom', \
        rotation=-20,fontsize=ccfs)
    caption += extract_short_caption(fname)


#########################################################################
### STIS Bar5 detection limit

if cfg['STIS']:
    fname = datapath+'HST_STIS.txt'
    a_STIS = ascii.read(fname)
    ax1.plot(a_STIS['Rho(as)'],a_STIS['KLIP_Contr'],color=c_bbvis,\
        linewidth=lw1,label='')
    ax1.text(0.2,5*10**-5,'HST STIS',color=c_bbvis,horizontalalignment='left',rotation=-40,fontsize=ccfs)
    caption += extract_short_caption(fname)


#########################################################################
### ACS detection limit

if cfg['ACS']:
    fname = datapath+'HST_ACS.txt'
    a_ACS = ascii.read(fname)
    ax1.plot(a_ACS['Rho(as)'],a_ACS['F606W_contr'],color=c_v,linewidth=lw1,label='')
    ax1.text(4,6*10**-9,'HST ACS',color=c_v,horizontalalignment='right',rotation=-35,fontsize=ccfs)
    caption += extract_short_caption(fname)


#########################################################################
### MagAO detection limit

if cfg['MagAO']:
    fname = datapath+'magao_ip_alphacen_5sigma.txt'
    a_MagAO_ip = ascii.read(fname)
    a_MagAO_ip['ip_Contrast'] = a_MagAO_ip['ip_contr_60min']
    ax1.plot(a_MagAO_ip['Rho(as)'], a_MagAO_ip['ip_Contrast'], \
        color=c_band4,linewidth=lw1,label='')
    ax1.plot([ a_MagAO_ip['Rho(as)'][-1], 1.7], \
        [0.9*a_MagAO_ip['ip_Contrast'][-1], 2E-8], 'k', linewidth=0.5)
    ax1.text(1.7, 2E-8,'Magellan VisAO',color=c_band4, horizontalalignment='left', \
        verticalalignment='top',rotation=-35,fontsize=ccfs)

    fname = datapath+'magao_Ys_betapic_5sigma.txt'
    a_MagAO_ys = ascii.read(fname)
    a_MagAO_ys['Ys_Contrast'] = a_MagAO_ys['Ys_contr_60min']
    ax1.plot(a_MagAO_ys['Rho(as)'], a_MagAO_ys['Ys_Contrast'], \
        color=c_yjh,linewidth=lw1,label='')
    ax1.text(1.7, 3.5E-7,'Magellan VisAO',color=c_yjh, horizontalalignment='left', \
        verticalalignment='bottom',rotation=-12,fontsize=ccfs)



#########################################################################
### SPHERE detection limit

if cfg['SPHERE']:
    fname = datapath+'SPHERE_Vigan.txt'
    a_SPHERE = ascii.read(fname)
    a_SPHERE['Contrast'] = 10**(-0.4*a_SPHERE['delta'])

    # manually split into IFS and IRDIS, at 0.7", as per documentation.
    idx_yjh = a_SPHERE['Rho(as)'] <= 0.7 # IFS YJH
    idx_k12 = a_SPHERE['Rho(as)'] >= 0.7  # IRDIS K1-K2
    ax1.plot(a_SPHERE['Rho(as)'][idx_yjh], a_SPHERE['Contrast'][idx_yjh], \
        color=c_yjh, linewidth=lw1, label='')
    ax1.plot(a_SPHERE['Rho(as)'][idx_k12], a_SPHERE['Contrast'][idx_k12], \
        color=c_k, linewidth=lw1, label='')
    ax1.text(0.19, 2.5E-6, 'VLT SPHERE', color=c_k, horizontalalignment='right', \
        verticalalignment='top', fontsize=ccfs)
    ax1.text(0.14, 1.3*10**-6, 'IFS /', color=c_yjh, horizontalalignment='right', \
        verticalalignment='top', fontsize=ccfs)
    ax1.text(0.14, 1.3*10**-6, ' IRDIS', color=c_k, horizontalalignment='left', \
        verticalalignment='top', fontsize=ccfs)
    caption += extract_short_caption(fname)


#########################################################################
### GPI H-band

if cfg['GPI']:
    fname = datapath+'GPI_Sirius_Ltype.txt'
    a_GPI = ascii.read(fname)
    ax1.plot(a_GPI['Rho(as)'],a_GPI['H_contr_60min_Ltype'],color=c_h,linewidth=lw1,label='')
    ax1.text(0.15,1.1*10**-5,'Gemini GPI',color=c_h,horizontalalignment='left',rotation=-24,fontsize=ccfs)
    caption += extract_short_caption(fname)


#########################################################################
### WFIRST


## predictions

if cfg['pred_img']:
    fname = datapath+'WFIRST_pred_imaging.txt'
    dat = ascii.read(fname)
    dat['lambda'].unit = u.nm
    dat['Rho(as)'] = dat['l/D'] * (dat['lambda'] / d_tel).decompose()*206265
    ax1.plot(dat['Rho(as)'], dat['contr_snr5'], color=c_v, linewidth=lw2, label='')
    ax1.text(dat['Rho(as)'][0], 1.5*dat['contr_snr5'][0], 'WFIRST  \nCGI pred.  ', color='darkblue',\
        horizontalalignment='right', verticalalignment='bottom', weight='bold', fontsize=ccfs+1)
    if cfg['exp_t']:
        if not cfg['req_img']:
            ax1.text(dat['Rho(as)'][0], 0.9*dat['contr_snr5'][-1], \
            'img, %ghr'%(dat['t_int_hr'][0]), color=c_v, weight='bold',\
            horizontalalignment='left', verticalalignment='top', fontsize=ccfs+1)
        else:
            ax1.text(dat['Rho(as)'][0], 0.9*dat['contr_snr5'][-1], \
            '%ghr'%(dat['t_int_hr'][0],), color=c_v, weight='bold',\
            horizontalalignment='left', verticalalignment='top', fontsize=ccfs+1)
    else:
        if not cfg['req_img']:
            ax1.text(dat['Rho(as)'][0], 0.9*dat['contr_snr5'][-1], 'img', color=c_v, weight='bold',\
                horizontalalignment='left', verticalalignment='top', fontsize=ccfs+1)
    caption += extract_short_caption(fname)


if cfg['pred_spec']:
    fname = datapath+'WFIRST_pred_spec.txt'
    dat = ascii.read(fname)
    dat['lambda'].unit = u.nm
    dat['Rho(as)'] = dat['l/D'] * (dat['lambda'] / d_tel).decompose()*206265
    dat['contr_snr5'] = dat['contr_snr10'] / 2.
    ax1.plot(dat['Rho(as)'], dat['contr_snr5'], color=c_band3, linewidth=lw2, label='')
    if not (cfg['pred_img'] or cfg['pred_disk']):
        ax1.text(dat['Rho(as)'][-1], dat['contr_snr5'][-1], ' WFIRST\n CGI pred.', color=c_band3,\
            horizontalalignment='left', verticalalignment='center', weight='bold', fontsize=ccfs+1)
    if cfg['exp_t']:
        if not cfg['req_spec']:
            ax1.text(dat['Rho(as)'][0], 1.1*dat['contr_snr5'][0], \
            'spec, %ghr'%(dat['t_int_hr'][0]), color=c_band3, weight='bold',\
            horizontalalignment='left', verticalalignment='bottom', fontsize=ccfs+1, zorder=6)
        else:
            ax1.text(dat['Rho(as)'][0], 1.1*dat['contr_snr5'][0], \
            '%ghr'%(dat['t_int_hr'][0]), color=c_band3, weight='bold',\
            horizontalalignment='left', verticalalignment='bottom', fontsize=ccfs+1, zorder=6)
    else:
        if not cfg['req_spec']:
            ax1.text(dat['Rho(as)'][0], 1.1*dat['contr_snr5'][0], \
            'spec', color=c_band3, weight='bold',\
            horizontalalignment='left', verticalalignment='bottom', fontsize=ccfs+1, zorder=6)
    caption += extract_short_caption(fname)



if cfg['pred_disk']:
    fname = datapath+'WFIRST_pred_disk.txt'
    dat = ascii.read(fname)
    ax1.plot(dat['Rho(as)'], dat['Band4_contr_snr5'], color=c_band4, linewidth=lw2, label='')
    if not cfg['pred_img']:
        ax1.text(1.3, 2E-9, 'WFIRST\nCGI pred.', color='darkred', horizontalalignment='left',\
            verticalalignment='center', fontsize=ccfs+1, weight='bold')
    if not cfg['req_disk']:
        ax1.text(dat['Rho(as)'][-1], 1.1*dat['Band4_contr_snr5'][-1], ' img, 100hr', color=c_band4,\
            horizontalalignment='left', verticalalignment='center', fontsize=ccfs+1, weight='bold')
    caption += extract_short_caption(fname)



## BTRs

if cfg['req_img']:
    fname = datapath+'WFIRST_req_imaging.txt'
    dat = ascii.read(fname)
    ax1.plot(dat['Rho(as)'], dat['Band1_contr_snr5'], color=c_v, linewidth=lw2, label='')
    ax1.text(dat['Rho(as)'][0], 1.1*dat['Band1_contr_snr5'][0], 'img', color=c_v,\
        horizontalalignment='left', verticalalignment='bottom', fontsize=ccfs+1)
    ax1.text(dat['Rho(as)'][0], dat['Band1_contr_snr5'][0], 'WFIRST \nCGI req. ', color='darkblue',\
        horizontalalignment='right', verticalalignment='center', weight='bold', fontsize=ccfs+1)
    caption += extract_short_caption(fname)

if cfg['req_disk']:
    fname = datapath+'WFIRST_req_disk.txt'
    dat = ascii.read(fname)
    ax1.plot(dat['Rho(as)'], dat['Band4_pt_contr_snr5'], color=c_band4, linewidth=lw2, label='')
    ax1.text(dat['Rho(as)'][-3], 1.1*dat['Band4_pt_contr_snr5'][-3], 'disk ', color=c_band4,\
        horizontalalignment='right', verticalalignment='bottom', fontsize=ccfs+1)
    if not (cfg['req_img'] or cfg['req_spec']):
        ax1.text(dat['Rho(as)'][0], dat['Band4_pt_contr_snr5'][0], ' WFIRST \n CGI req. ', color='darkblue',\
            horizontalalignment='right', verticalalignment='center', weight='bold', fontsize=ccfs+1)
    caption += extract_short_caption(fname)

if cfg['req_spec']:
    import matplotlib.patheffects as path_effects
    fname = datapath+'WFIRST_req_spec.txt'
    dat = ascii.read(fname)
    if cfg['req_disk']:  # draw a shadow under the line to make it easier to see if overlapping other req line
        p = [path_effects.SimpleLineShadow(offset=(1, -1)), path_effects.Normal()]
    else:
        p = [path_effects.Normal()]
    ax1.plot(dat['Rho(as)'], dat['Band3_contr_snr5'], color=c_band3, linewidth=lw2-0.5, \
        label='', path_effects=p)
    ax1.text(dat['Rho(as)'][0], 1.1*dat['Band3_contr_snr5'][0], 'spec', color=c_band3,\
        horizontalalignment='left', verticalalignment='bottom', fontsize=ccfs+1)
    if not cfg['req_img']:
        ax1.text(dat['Rho(as)'][0], dat['Band3_contr_snr5'][0], ' WFIRST  \n CGI req.  ', color='darkblue',\
            horizontalalignment='right', verticalalignment='center', weight='bold', fontsize=ccfs+1)
    caption += extract_short_caption(fname)



#########################################################################
###### -------------------- planets -------------------------  ##########
#########################################################################

#########################################################################
### Self luminous directly imaged planets

if cfg['DI_H'] or cfg['DI_B1_pred'] or cfg['DI_B3_pred']:
    fname = datapath+'DIplanets.txt'
    a_DI = ascii.read(fname)
    a_DI['547m_contr'] = 10**(a_DI['547m_delta']/-2.5)
    a_DI['763m_contr'] = 10**(a_DI['763m_delta']/-2.5)
    a_DI['H_contr'] = 10**(a_DI['H_delta']/-2.5)
    alpha_di = 0.7
    caption += extract_short_caption(fname)

if cfg['DI_H']:
    ax1.scatter(a_DI['Rho(as)'],a_DI['H_contr'],color=c_h, edgecolor='k', \
        alpha=alpha_di, marker='s', s=cfg['di_markersize']-15, zorder=2,\
        label='directly imaged, 1.6$\mathdefault{\mu} $m observed')

if cfg['DI_B3_pred']:
    ax1.scatter(a_DI['Rho(as)'],a_DI['763m_contr'],color=c_band3, edgecolor='k', \
        marker='d', alpha=alpha_di, s=cfg['di_markersize'], zorder=2, \
        label='directly imaged, 750nm predicted')
    if not cfg['DI_B1_pred']:
        for ct, rho in enumerate(a_DI['Rho(as)']):
            ax1.plot([rho,rho], [a_DI[ct]['763m_contr'], a_DI[ct]['H_contr']], \
            color='lightgray', linewidth=1, linestyle=':', zorder=1)

if cfg['DI_B1_pred']:
    for ct, rho in enumerate(a_DI['Rho(as)']):
        ax1.plot([rho,rho], [a_DI[ct]['547m_contr'], a_DI[ct]['H_contr']], \
            color='lightgray', linewidth=1, linestyle=':', zorder=1)
    ax1.scatter(a_DI['Rho(as)'],a_DI['547m_contr'],color=c_v, edgecolor='k', \
        marker='o', alpha=alpha_di, s=cfg['di_markersize'], zorder=2, \
        label='directly imaged, 550nm predicted')


#########################################################################
### Specific planetary systems

if cfg['Prox_Cen']:
    albedo = 0.35
    sma = 0.05*u.au
    flux_ratio = rlp.calc_lambert_flux_ratio(sma=sma, rp=1.3**(1./3)*u.earthRad,\
        orb_ang=0*u.degree, albedo=albedo, inclin=0*u.degree)
    rho = (sma/1.3*u.pc).value
    ax1.scatter(rho,flux_ratio,marker='^', color=c_pl)
    ax1.text(rho,flux_ratio,'  Proxima Cen b',color=c_pl,\
        horizontalalignment='left',verticalalignment='center',fontsize=ccfs)
    caption += '-- Proxima Cen b. At quadrature, albedo = ' + str(albedo) +\
        ', radius = (M/Me)^(1/3) * Re, circular orbits.\n\n'


if cfg['Tau_Ceti']:
    # Tau Ceti
    tc_dist = 3.65*u.pc
    albedo = 0.35
    # make a separate upper x axis for physical separation of this system
    ax3 = ax1.twiny()
    ax3.set_ylim(ylim)
    ax3.set_xlim(xlim * tc_dist.value)
    ax3.set_xscale('log')
    if cfg['Tau_Ceti_axis']:
        ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax3.set_xlabel('Semi-major axis for Tau Ceti (d=%.3g pc) [au] '%tc_dist.value)
    else:
        ax3.set_xticks([])
        ax3.xaxis.set_ticks_position('none')

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
        horizontalalignment='right',verticalalignment='top',fontsize=ccfs)

    caption += '-- Tau Ceti e&f. At quadrature, albedo = ' + str(albedo) +\
        ', radius = (M/Me)^(1/3) * Re, circular orbits.\n\n'


if cfg['solar_system']:
    # Earth & Jupiter
    #earthRatio = rlp.calc_lambert_flux_ratio(sma=1.*u.au, rp=1.*u.earthRad,\
    #    orb_ang=0*u.degree,albedo=0.367, inclin=0*u.degree)
    earthRatio = 1.0E-10 # use the "standard" value
    jupiterRatio = rlp.calc_lambert_flux_ratio(sma=5.*u.au, rp=1.*u.jupiterRad,\
        orb_ang=0*u.degree,albedo=0.52, inclin=0*u.degree)
    caption += '-- Earth (==1E-10) & Jupiter (albedo=0.52) at quadrature as seen from 10 pc. '+\
                '(Jupiter albedo: Traub & Oppenheimer, '+\
                'Direct Imaging chapter of Seager Exoplanets textbook, Table 3)\n\n'

    ax1.scatter(0.1,earthRatio,marker='$\\bigoplus$',color=c_pl, s=cfg['rv_markersize'], zorder=5)
    ax1.text(0.1,earthRatio,'  Earth at 10pc',color=c_pl,\
        horizontalalignment='left',verticalalignment='center',fontsize=ccfs)

    ax1.scatter(0.5,jupiterRatio,marker='v',color=c_pl,  edgecolor='k', s=cfg['rv_markersize'], zorder=5)
    ax1.text(0.5,jupiterRatio,' Jupiter at 10pc',color=c_pl,\
    horizontalalignment='left',verticalalignment='top',fontsize=ccfs)

#    ax1.text(xlim[1], ylim[0]*1.1, 'Solar System as seen from 10pc. ',\
#    color='k',horizontalalignment='right', verticalalignment='bottom',fontsize=ccfs-1)


#########################################################################
###Add RV planets

if cfg['RV_pred']:
    fname = datapath+'reflected_light_table.txt'
    tmp = ascii.read(fname)
    idx_rv = tmp['pl_discmethod'] == "Radial Velocity"
    ax1.scatter(tmp[idx_rv]['sma_arcsec'], tmp[idx_rv]['Fp/F*_quad'], s=cfg['rv_markersize'],\
        color='dimgray', edgecolor='k', marker='^', label='RV, reflected light, predicted', zorder=5)
    caption += extract_short_caption(fname)


#########################################################################
###Plot axes, tick mark ajdusting, legend, etc.

if cfg['is_draft']:
    ax1.text(xlim[1], ylim[0]*2, "DRAFT  "+str(date.today()) + ' ', \
        horizontalalignment='right',verticalalignment='bottom',color='red',weight='bold')

first_legend = ax1.legend(fontsize=cfg['legend_font_size'], loc='upper right', \
    title='Known Exoplanets')
first_legend.get_title().set_fontsize(8)

ax1.grid(b=True, which='major', color='tan', linestyle='-', alpha=0.1)

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xticks([0.05, 0.1, 0.5, 1, 5])
ax1.set_ylim(ylim)
ax1.set_xlim(xlim)

# write x axis in scalar notation instead of powers of 10
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1g'))

if cfg['color_by_lambda'].lower() != 'none':
    second_legend = ax2.legend(loc='upper left', fontsize=cfg['legend_font_size'], title='Wavelength ($\lambda_0$)')
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

try:
    os.stat('./output')
except:
    os.mkdir('./output')


with open('./output/auto_caption_'+cfgname+'.txt','w') as f:
    f.write(caption)

if cfg['save_pdf']:
    plt.savefig('./output/flux_ratio_'+cfgname+'.pdf')
if cfg['save_jpg']:
    plt.savefig('./output/flux_ratio_'+cfgname+'.jpg', dpi=cfg['jpg_dpi'])
if cfg['save_png']:
    plt.savefig('./output/flux_ratio_'+cfgname+'.png', dpi=cfg['png_dpi'])
