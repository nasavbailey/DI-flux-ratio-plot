import matplotlib.pyplot as plt
import numpy as np

def calc_alpha(i_dgr, orb_ang_dgr, out_dgr=True):
    #Calculate the phase angle (star/planet/observer angle) for a
    # given system inclination and orbit angle. Assumes a *circular* orbit
    #
    #Inputs:
    #i_dgr = inclination of orbit, in degrees [single value]
    #orb_ang_dgr = angle in the orbit, relative to the ascending node [single value or vector]
    #            = true anomaly + argument of periastron
    #out_dgr = [True]/False -- output alpha in degrees?
    #
    #References
    # http://iopscience.iop.org/article/10.1088/0004-637X/747/1/25; Madhusudhan & Burrows 2012
    dgr2rad = np.pi/180.0
    alpha = np.arccos( np.sin(orb_ang_dgr*dgr2rad) * np.sin(i_dgr*dgr2rad) )
    if out_dgr:
        return alpha / dgr2rad
    else:
        return alpha

def calc_lambert(alpha_rad):
    if max(np.abs(alpha_rad)) > 2*np.pi:
        raise Exception('alpha must be in radians')
    phi = ( np.sin(alpha_rad) + (np.pi - alpha_rad) * np.cos(alpha_rad) ) / np.pi
    return phi

def calc_xy_r_int(i_dgr):
    a_dgr = np.linspace(0,360,50)
    a_rad = a_dgr * np.pi/180
    alpha = calc_alpha(i_dgr,a_dgr,out_dgr=False)
    inten = calc_lambert(alpha)
    x = np.cos(a_rad)
    y = np.sin(a_rad) * np.cos(i_dgr*np.pi/180)
    r = np.sqrt(x**2 + y**2)
    return (x,y,r,inten)

def max_flux:
    return
#https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets&select=pl_hostname,st_dist,st_vj,st_ic,pl_orbsmax,pl_orbincl,pl_bmassj&order=pl_hostname&format=ascii
#https://exoplanetarchive.ipac.caltech.edu/docs/API_exoplanet_columns.html#defaultcol
#http://mips.as.arizona.edu/~cnaw/sun.html

inclination = 60
(x,y,r,inten) = calc_xy_r_int(inclination)


ang = np.linspace(90,270,50)
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)

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

#ax3.plot(alpha*180/np.pi, calc_lambert(alpha))
#ax3.set_xlabel('alpha')
#ax3.set_ylabel('flux relative to full')

ax3.plot(r, inten)
ax3.set_ylim([0,1])
ax3.set_xlabel('normalized radius')
ax3.set_ylabel('normalized flux'    )

ax4.scatter(x,y, c=inten, vmin=0, vmax=1 )
ax4.set_xlabel('x position')
ax4.set_xlabel('y position')
ax4.axis('equal')

plt.tight_layout()
plt.savefig('rv.pdf')

print 'done'
