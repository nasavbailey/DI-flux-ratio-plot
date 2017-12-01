import matplotlib.pyplot as plt
import numpy as np

def calc_alpha(i_dgr, orb_ang_dgr, out_dgr=True):
    #Calculate the phase angle for a given system inclination and orbit angle
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


ang = np.linspace(90,270,50)
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
for inc in np.linspace(0,90,5):
    alpha = calc_alpha(inc,ang,out_dgr=False)
    ax1.plot(ang, alpha*180/np.pi)
    ax2.plot(ang, calc_lambert(alpha), label=inc)

ax3.plot(alpha*180/np.pi, calc_lambert(alpha))


ax1.set_xlabel('angle from ascending node [dgr]')
ax1.set_ylabel('alpha')
ax1.legend()

ax2.axhline(0.25,color='gray',linestyle='-.')
ax2.axhline(1/np.pi,color='lightgray',linestyle=":")
ax2.set_xlabel('angle from ascending node [dgr]')
ax2.set_ylabel('flux relative to full')
ax2.legend(title='inclination', fontsize=6)

ax3.axhline(0.25,color='gray',linestyle='-.')
ax3.axhline(1/np.pi,color='lightgray',linestyle=":")

plt.tight_layout()
plt.savefig('rv.pdf')

print 'done'
