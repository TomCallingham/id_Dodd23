from .Potential_H99 import H99_potential

import numpy as np
def dynamics_calc_H99(xyz:np.ndarray)->np.ndarray:
    '''
    In: [N,6] np array of [x,y,z,vx,vy,vz]
    Out: [N,3] np array of [E,Lz,Lperp]
    '''
    pos, vel = xyz[:,:3], xyz[:,3:]
    K = 0.5 * np.sum(vel**2, axis=1)
    U = H99_potential(pos)
    En = U+K
    Lvec = np.cross(vel, pos)
    Lz = Lvec[:, 2]
    Lperp = np.sqrt((Lvec[:,:2]**2).sum(axis=-1))
    vec_ELzLp = np.column_stack((En,Lz,Lperp))
    return vec_ELzLp
