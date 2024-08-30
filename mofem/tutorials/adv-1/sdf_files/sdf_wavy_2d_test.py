import numpy as np
A = 0.0002
wave_len = 2.
w = 2. * np.pi / wave_len
ind = 0.00159624

class WavySurface:
    def sDF(x, y,z, t):
        sdf =  (A * (1. - np.cos(2. * np.pi * x / wave_len)) -  y) -ind * t
        return sdf

    def gradSDF(x, y, z, t):
        df_dx = A * (2. * np.pi /wave_len) * np.sin(2. * np.pi * x / wave_len)
        df_dy = np.full_like(df_dx, -1.0)
        df_dz = np.zeros_like(df_dx)
        return np.hstack([df_dx.reshape((-1, 1)), df_dy.reshape((-1, 1)), df_dz.reshape((-1, 1))])

    def hessSDF(x, y, z, t):
        fxx = A * ((2. * np.pi / wave_len)**2) * np.cos(2. * np.pi * x / wave_len)
        fxy = np.zeros_like(fxx)
        fyx = np.zeros_like(fxx)
        fyy = np.zeros_like(fxx)
        fzx = np.zeros_like(fxx)
        fzy = np.zeros_like(fxx)
        fzz = np.zeros_like(fxx)
        return np.hstack([fxx.reshape((-1,1)), fxy.reshape((-1,1)), fzx.reshape((-1,1)),fyy.reshape((-1,1)), fzy.reshape((-1,1)), fzz.reshape((-1,1))])

def sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return WavySurface.sDF(x, y, z, t)

def grad_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return WavySurface.gradSDF(x, y,z, t)

def hess_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return WavySurface.hessSDF(x, y,z, t)