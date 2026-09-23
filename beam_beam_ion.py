import numpy as np
from scipy.special import wofz


def gaussian_be_field(separation, capcov):
    """Bassetti-Erskine K=(Kx,Ky) in 1/mm for 2D Gaussian capcov [mm^2]."""
    d = np.asarray(separation, dtype=float)
    eig, rotation = np.linalg.eigh(np.asarray(capcov, dtype=float))
    if eig[0] <= 0:
        raise ValueError('Convolved covariance must be positive definite')
    eig = eig[::-1]
    rotation = rotation[:, ::-1]
    u, v = rotation.T @ d
    sx, sy = np.sqrt(eig)
    r2 = u*u + v*v
    if r2 == 0:
        return np.zeros(2)
    if abs(sx-sy)/sx < 1e-6:
        field = -np.expm1(-r2/(2*sx*sx))/r2 * np.array([u, v])
    elif r2 < 1e-10 * sy*sy:
        field = np.array([u/sx, v/sy])/(sx+sy)
    else:
        a, b = abs(u), abs(v)
        scale = np.sqrt(2*(sx*sx-sy*sy))
        exponent = -a*a/(2*sx*sx) - b*b/(2*sy*sy)
        w = -1j*np.sqrt(np.pi)/scale * (
            wofz((a+1j*b)/scale) - np.exp(exponent)*
            wofz((a*sy/sx+1j*b*sx/sy)/scale))
        field = np.array([np.sign(u)*w.real, -np.sign(v)*w.imag])
    return rotation @ field


def ion_bb_displacements(separation, N1_ions, N2_ions, Z, A, gamma, capcov,
                         beta_star, tunes, doros_arm, proton_radius):
    """Return [B1,B2] x [X,Y] displacements in mm for IP and DOROS."""
    K = gaussian_be_field(separation, capcov)
    # Javier/Marcus sign convention, d = set(B1)-set(B2).
    theta1 = -2*N2_ions*(Z*Z/A)*proton_radius/gamma * K
    theta2 = +2*N1_ions*(Z*Z/A)*proton_radius/gamma * K
    arc_factor = beta_star/(2*np.tan(np.pi*np.asarray(tunes)))
    doros_factor = arc_factor + doros_arm/2
    angles = np.stack([theta1, theta2])
    return angles*arc_factor, angles*doros_factor
