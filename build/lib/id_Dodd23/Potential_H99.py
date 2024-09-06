
import numpy as np

################################################################################
# ORIGINAL FILE BY JOVAN VELJANOSKI
################################################################################

G = 4.301*10**-6   # fixed GMsun [kpc km^2/s^2]

M200 = 1e12
rs = 21.5
Rvir = 258

disc_a = 6.5
disc_b = 0.26
Mdisc = 9.3*10**+10

bulge_c = 0.7
Mbulge = 3.0*10**+10

def H99_potential(pos):
    #pos galactic position, kpc, Nx3
    phi_halo  = potential_halo(pos) #NFW
    phi_disc  = potential_disc(pos) #Nig?
    phi_bulge = potential_bulge(pos) #standard
    phi_total = phi_disc + phi_bulge + phi_halo
    return phi_total
################################################################################

def potential_halo(pos):
    c = Rvir/rs
    phi_0 = G*M200/Rvir / (np.log(1+c)-c/(1+c))*c
    r = np.linalg.norm(pos,axis=1)
    phi_halo = - phi_0 * rs/r * np.log(1.0 + r/rs)
    return phi_halo



def potential_disc(pos):
    '''
    Calculates the potential due to the disc
    '''
    # parameters
    GMd = G * Mdisc

    x, y, z = pos[:,0], pos[:,1], pos[:,2]

    sqd = np.sqrt(z**2.0 + disc_b**2.0)
    # square root of the density, probably
    sqden1 = np.sqrt(x**2. + y**2.0 + (disc_a+sqd)**2.0)
    # the potential of the disc
    phi_d = -GMd/sqden1

    return phi_d

################################################################################


def potential_bulge(pos):
    '''
    Calculates the potential contribution due to the bulge
    '''
    # parameters
    GMb = G * Mbulge

    # radial distance
    r = np.linalg.norm(pos,axis=1)

    # the potential due to the bulge
    phi_b = -GMb/(r+bulge_c)

    return phi_b

################################################################################



