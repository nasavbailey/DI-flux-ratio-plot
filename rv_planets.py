#import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.io import ascii


def calc_alpha(i_dgr, orb_ang_dgr, eccen=0, arg_peri_dgr=0, out_dgr=True):
    # Calculate the phase angle (star/planet/observer angle) for a
    #  given system inclination and orbit angle.
    # Assumes a *circular* orbit
    #
    # Inputs:
    # i_dgr = inclination of orbit, in degrees [single value]
    # orb_ang_dgr = angle in the orbit, relative to the ascending node [single value or vector]
    #            = true anomaly + argument of periastron
    # out_dgr = [True]/False -- output alpha in degrees?
    #
    # Future work:
    # eccen = orbit eccentricity. NOT implemented yet.
    # arg_peri_dgr = argument of periastron in degrees. NOT implemented yet.
    #
    # References
    # http://iopscience.iop.org/article/10.1088/0004-637X/747/1/25; Madhusudhan & Burrows 2012

    if eccen !=0 or arg_peri_dgr != 0:
        raise Exception('Non-circular orbits not yet implemented')

    alpha = np.arccos( np.sin(orb_ang_dgr*u.degree) * np.sin(i_dgr*u.degree) )
    if out_dgr:
        return alpha / dgr2rad
    else:
        return alpha


def calc_lambert(alpha_dgr):
    # Evaluate the phase law of a Lambertian planet at a given phase angle
    #
    # Inputs:
    # Phase angle in degrees
    #
    # Outputs:
    # Phi - a value between 0 and 1, where 1 is the value at alpha=0, and 0 is no flux
    #
    # References:
    # Traub & Oppenheimer: Direct Imaging of Exoplanets, equation 15.
    #  (pg 116 of Seager Exoplanets textbook)

    alpha_rad = alpha_dgr * np.pi/180.
    phi = ( np.sin(alpha_rad) + (np.pi - alpha_rad) * np.cos(alpha_rad) ) / np.pi
    return phi


def calc_xy_r_int(i_dgr):
    # Calculate the xy positions, projected separation, & lambertian brightness
    #  of a planet on a CIRCULAR, inclined orbit
    #
    # Inputs:
    # i_dgr = inclination in degrees
    #
    # Outputs:
    # list of vectors (x,y,r,inten)

    ang_dgr = np.linspace(0,360,70)
    ang_rad = a_dgr * np.pi/180
    alpha_dgr = calc_alpha(i_dgr,ang_dgr,out_dgr=True)
    inten = calc_lambert(alpha_dgr)
    x = np.cos(ang_rad)
    y = np.sin(ang_rad) * np.cos(i_dgr*np.pi/180)
    r = np.sqrt(x**2 + y**2)
    return (x,y,r,inten)


def lambertian_phase_curve(rho_au, rp=1.*u.jupiterRad, albedo=0.5, \
        incl_dgr=90, orb_ang_dgr=0):
    # Calculate flux ratio of a Lambertian planet on a CIRCULAR orbit

    phi = calc_lambert(calc_alpha(incl_dgr, orb_ang_dgr, out_dgr=False))

    rad_ratio = ((rp / rho_au)**2).decompose() # pl_orbsmax already has associated unit
    f_ratio = albedo * phi * rad_ratio
    return f_ratio


def read_and_filter_exo_archive(fname='exo_archive_query.txt', \
        rho_min=0.14, rho_max=10., st_v_min=8., mp_min=0.25):
    dat = ascii.read(fname)

    # remove rows with missing important values
    filter_cols = ['st_dist','st_vj', 'pl_orbsmax','pl_bmassj']
    for f in filter_cols:
        dat = dat[~dat[f].mask]

    # replace spaces with _ in star names
    dat['pl_hostname'] = [x.replace(' ','_') for x in dat['pl_hostname']]

    # keep only planets with a maximum elongation between X1 and X2
    dat['a_as'] = dat['pl_orbsmax'] / dat['st_dist']
    dat.sort('a_as')
    dat['a_as'].format = '%.2f'
    dat['a_as'].unit = u.arcsec
    dat = dat[dat['a_as'] >= rho_min]
    dat = dat[dat['a_as'] < rho_max]

    #keep only stars brighter than st_v_min in V
    dat = dat[dat['st_vj'] < st_v_min]

    #keep only planets larger than mp_min*Mj
    dat = dat[dat['pl_bmassj'] > mp_min]

    return dat


def demo_phase():
    plt.figure(figsize=[8,4])
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)

    ang = np.linspace(90,270,50)
    for inc in np.linspace(0,90,5):
        alpha = calc_alpha(inc,ang,out_dgr=False)
        ax1.plot(ang, alpha*180/np.pi)
        ax2.plot(ang, calc_lambert(alpha), label=inc)

    ax1.set_xlabel('angle from ascending node [dgr]')
    ax1.set_ylabel('alpha')
    ax1.legend()

    ax2.set_xlabel('angle from ascending node [dgr]')
    ax2.set_ylabel('flux relative to full')
    ax2.legend(title='inclination', fontsize=6)

    (x,y,r,inten) = calc_xy_r_int(60)
    ax3.plot(r,inten,label='i=60')
    #ax1.set_yscale('log')
    plt.legend()

    ax4.scatter(x,y,c=inten, vmin=0, vmax=1,label='i=60')
    ax4.axis('equal')
    plt.legend()

    plt.tight_layout()
    plt.savefig('demo_phase_curve.pdf')
