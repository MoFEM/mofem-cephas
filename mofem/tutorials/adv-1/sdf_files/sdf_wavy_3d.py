import numpy as np

A = 0.0002
wave_len = 2.0
w = 2. * np.pi / wave_len
ind = 0.0021745

class WavySurface:
    def sDF(x, y,z, t):
        sdf =  A * (1. - np.cos(2 * np.pi * x / wave_len) * np.cos(2. * np.pi * y / wave_len)) - z - ind * t
        return sdf

    def gradSDF(x, y, z, t):
        df_dx = A * (2. * np.pi / wave_len) * np.sin(2. * np.pi * x / wave_len) * np.cos(2. * np.pi * y / wave_len)
        df_dy = A * (2. * np.pi / wave_len) * np.cos(2 * np.pi * x / wave_len) * np.sin(2. * np.pi * y / wave_len)
        df_dz = np.full_like(df_dx, -1.0)
        return np.hstack([df_dx.reshape((-1, 1)), df_dy.reshape((-1, 1)), df_dz.reshape((-1, 1))])

    def hessSDF(x, y, z, t):
        # symetric hess mat pattern: xx, xy, xz, yy, yz, zz
        fxx = A * ((2. * np.pi / wave_len)**2) * np.cos(2. * np.pi * x / wave_len) * np.cos(2. * np.pi * y / wave_len)
        fxy = - A * ((2 * np.pi / wave_len)**2) * np.sin(2 * np.pi * x / wave_len) * np.sin(2 * np.pi * y / wave_len)
        fzx = np.zeros_like(fxx)
        fyy =  A * ((2. * np.pi / wave_len)**2) * np.cos(2. * np.pi * x / wave_len) * np.cos(2. * np.pi * y / wave_len)
        fyz = np.zeros_like(fyy)
        fzz = np.zeros_like(fxx)
        return np.hstack([fxx.reshape((-1,1)), fxy.reshape((-1,1)), fzx.reshape((-1,1)),fyy.reshape((-1,1)), fyz.reshape((-1,1)), fzz.reshape((-1,1))])

def sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return WavySurface.sDF(x, y, z, t)

def grad_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return WavySurface.gradSDF(x, y,z, t)

def hess_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return WavySurface.hessSDF(x, y,z, t)