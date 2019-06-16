#import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.io import ascii
from astropy.table import QTable

def unit_check_convert(var, myunit):
    # Check variable units.
    # Add unit if none. If mismatch: convert if possible; error if incompatible.
    #
    # Inputs:
    # var = variable to check. Single value or array.
    # myunit = desired unit in astropy unit format (eg: u.degree)
    #
    # Output:
    # astropy quantity with specified unit. If any quantity in an input array can't be converted,
    #   will error even if other entries are OK.

    try:
        myunit.decompose()
    except:
        raise Exception('to_unit input must be an astropy unit, such as u.degree. You input: "%s" of %s'\
            % (str(myunit), str(type(myunit))))

    # tack on unit if var doesn't have units yet
    try:
        var.unit
    except:
        return var * myunit

    # check whether var units match the requested ones.
    if var.unit == myunit:
        return var
    else:
        print("...Converting [%s] to [%s]" % (str(var.unit), str(myunit)))
        return var.to(myunit)


def calc_alpha(i_dgr, orb_ang_dgr, eccen=0, arg_peri_dgr=0, out_dgr=True):
    # Calculate the phase angle (star/planet/observer angle) for a
    #  given system inclination and orbit angle.
    # Assumes a *circular* orbit
    #
    # Inputs:
    # i_dgr = inclination of orbit, in degrees [single value]
    # orb_ang = angle in the orbit, relative to the ascending node [single value or vector]
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

    i_dgr = unit_check_convert(i_dgr, u.degree)
    orb_ang_dgr = unit_check_convert(orb_ang_dgr, u.degree)

    alpha = np.arccos( np.sin(orb_ang_dgr) * np.sin(i_dgr) )
    if out_dgr:
        return alpha.to(u.degree)
    else:
        return alpha


def calc_lambert_phase_law(alpha):
    # Evaluate the phase law of a Lambertian planet at a given phase angle
    #
    # Inputs:
    # Phase angle in DEGREES
    #
    # Outputs:
    # Phi - a value between 0 and 1, where 1 is the value at alpha=0, and 0 is no flux
    #
    # References:
    # Traub & Oppenheimer: Direct Imaging of Exoplanets, equation 15.
    #  (pg 116 of Seager Exoplanets textbook)

    alpha = unit_check_convert(alpha, u.degree)

    return ( np.sin(alpha) + (np.pi - alpha.to(u.radian).value) * np.cos(alpha) ) / np.pi


def calc_lambert_flux_ratio(sma=3*u.au, rp=1.*u.jupiterRad, albedo=0.5, \
        inclin=90*u.degree, orb_ang=0*u.degree):
    # Calculate flux ratio of a Lambertian planet on a CIRCULAR, inclined orbit
    #
    # Inputs [defaults in brackets]:
    # rho     = orbital separation in AU (physical, NOT projected); [3au]
    # rp      = planet radius, in Jupiter radii; [1Rj]
    # albedo  = geometric albedo; [0.5]
    # incl    = inclination in degrees; [90dgr]
    # orb_ang = angle in the orbit, relative to the ascending node (single value or vector)
    #            == true anomaly + argument of periastron;  [0dgr = quadrature]
    # Outputs:
    # flux ratio

    sma = unit_check_convert(sma, u.au)
    rp = unit_check_convert(rp, u.jupiterRad)
    inclin = unit_check_convert(inclin, u.degree)
    orb_ang = unit_check_convert(orb_ang, u.degree)

    phi = calc_lambert_phase_law(calc_alpha(inclin, orb_ang))
    rad_ratio = ((rp / sma).decompose())**2 # pl_orbsmax already has associated unit
    f_ratio = albedo * phi * rad_ratio
    return f_ratio


def calc_xy_r_int(inclin):
    # Calculate the xy positions, projected separation, & lambertian brightness
    #  of a planet on a CIRCULAR, inclined orbit
    #
    # Inputs:
    # incl = inclination in degrees
    #
    # Outputs:
    # list of vectors (x,y,r,inten)

    inclin = unit_check_convert(inclin, u.degree)

    ang = np.linspace(0,360,70) * u.degree

    alpha = calc_alpha(inclin,ang,out_dgr=True)
    alpha = unit_check_convert(alpha, u.degree)

    inten = calc_lambert_phase_law(alpha)
    x = np.cos(ang)
    y = np.sin(ang) * np.cos(inclin)
    r = np.sqrt(x**2 + y**2)
    return (x,y,r,inten)


def read_and_filter_exo_archive(fname='data/exo_archive_query.txt', \
        rho_min=0.1*u.arcsec, rho_max=1.5*u.arcsec,\
        st_v_min=8.*u.mag, mp_min=0.25*u.jupiterMass):

    # make sure there are no conflicting units
    rho_min = unit_check_convert(rho_min, u.arcsec)
    rho_max = unit_check_convert(rho_max, u.arcsec)
    st_v_min = unit_check_convert(st_v_min, u.mag)
    mp_min = unit_check_convert(mp_min, u.jupiterMass)

    dat = ascii.read(fname)

    # remove rows with missing important values
    filter_cols = ['st_dist','st_optmag', 'pl_orbsmax','pl_bmassj']
    for f in filter_cols:
        dat = dat[~dat[f].mask]

    # manually assign b/c table units were (non-standard) "mags" not "mag"
    # Must be done before converting to QTable
    dat['st_optmag'].unit = u.mag

    dat = QTable(dat)

    # remove the metadata that had column descriptions, since it makes it a pain
    #  to write out data later if we change columns. Better to add column descriptions
    #  later by hand if needed...
    dat.meta = {}

    dat.rename_column('pl_orbsmax','sma_au')
    dat.rename_column('pl_orbincl','orb_incl')
    dat.rename_column('pl_orbeccen','eccen')
    dat.rename_column('pl_bmassj','pl_massj')

    # replace spaces with _ in star names
    dat['pl_hostname'] = [x.replace(' ','_') for x in dat['pl_hostname']]

    #keep only stars brighter than st_v_min in V
    dat = dat[dat['st_optmag'] < st_v_min]

    #keep only planets larger than mp_min*Mj
    dat = dat[dat['pl_massj'] > mp_min]

    # keep only planets with a maximum elongation between X1 and X2
    tmp = (dat['sma_au'] / dat['st_dist']).decompose() * u.radian
    dat['sma_arcsec'] = tmp.to(u.arcsec)
    dat['sma_arcsec'].info.format = '%.2f'
    dat = dat[dat['sma_arcsec'] >= rho_min]
    dat = dat[dat['sma_arcsec'] < rho_max]
    dat.sort('sma_arcsec')

    return dat



def create_wfirst_reflected_light_table(fname_in='data/exo_archive_query.txt', \
        fname_out = 'test.txt', \
        rho_min=0.12*u.arcsec, rho_max=1.45*u.arcsec,\
        st_v_min=7.*u.mag, mp_min=0.25*u.jupiterMass, rp=1.0*u.jupiterRad,\
        albedo=0.5, orb_ang=0*u.degree, inclin=90*u.degree):


    rho_min = unit_check_convert(rho_min, u.arcsec)
    rho_max = unit_check_convert(rho_max, u.arcsec)
    st_v_min = unit_check_convert(st_v_min, u.mag)
    mp_min = unit_check_convert(mp_min, u.jupiterMass)
    rp = unit_check_convert(rp, u.jupiterRad)
    inclin = unit_check_convert(inclin, u.degree)
    orb_ang = unit_check_convert(orb_ang, u.degree)

    # read the NASA exoplanet archive table & select appropriate systems
    dat = read_and_filter_exo_archive(fname=fname_in,rho_min=rho_min, rho_max=rho_max,\
            st_v_min=st_v_min, mp_min=mp_min)

    # Calculate the Lambertian flux ratio
    dat['Fp/F*_quad'] = calc_lambert_flux_ratio(sma=dat['sma_au'], \
        rp=rp,orb_ang=orb_ang,albedo=albedo, inclin=inclin)
    dat['Fp/F*_quad'].info.format = '%.1e'


    # create header comments / description
    comments = ('#Short caption: All planets from NASA exoplanet archive with a '+ \
        'semi-major axis of %.2f-%.2g arcsec, mass > %.2f Mjup, and host star V mag < %.2f. '+\
        'Lambertian flux ratio assumes: radius = %.1g Rjup, geometric albedo = %.2g, '+\
        'circular orbit, inclination = %.1f, and angle of %.1f degrees from the ascending node.\n' )% \
        (rho_min.value, rho_max.value, mp_min.value, st_v_min.value, rp.value, albedo, \
        inclin.to(u.degree).value, orb_ang.to(u.degree).value)
    comments += '#References:\n'
    comments += '#Lambertian phase curve from Traub & Oppenheimer: ' + \
    'Direct Imaging of Exoplanets, equation 15, pg 116 of Seager Exoplanets textbook.\n'
    #print comments

    #write out with astropy tables to format the body correctly
    t2 = dat['pl_hostname','pl_letter','sma_arcsec','Fp/F*_quad','st_optmag','st_spstr','pl_discmethod']
    ascii.write(t2,fname_out,\
            format='fixed_width', comment='#',delimiter=',',bookend=False,overwrite=True)

    # Then redo write out to get the header formatted the same as the other datafiles
    f = open(fname_out,'r')
    lines = f.readlines()
    f.close()
    f = open(fname_out,'w')
    f.write('#'+lines[0])
    f.write(comments)
    f.write(('').join(lines[1:]))
    f.close()

    print('File written to '+fname_out)

    return


def demo_phase():
    import matplotlib.pyplot as plt

    plt.figure(figsize=[6,6])
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)

    ang = np.linspace(90,270,50)*u.degree
    for inc in np.linspace(0,90,5)*u.degree:
        alpha = calc_alpha(inc,ang)
        ax1.plot(ang, alpha)
        ax2.plot(ang, calc_lambert_phase_law(alpha), label=inc)

    ax1.set_xlabel('angle from ascending node [dgr]')
    ax1.set_ylabel('alpha [dgr]')
    ax1.legend()

    ax2.set_xlabel('angle from ascending node [dgr]')
    ax2.set_ylabel('flux relative to full')
    ax2.legend(title='inclination', fontsize=6)

    (x,y,r,inten) = calc_xy_r_int(60)
    ax3.plot(r,inten,label='i=60')
    plt.legend()

    ax4.scatter(x,y,c=inten, vmin=0, vmax=1,label='i=60')
    ax4.axis('equal')
    plt.legend()

    plt.tight_layout()
    plt.savefig('demo_phase_curve.pdf')
    print('Saved plot called demo_phase_curve.pdf')
