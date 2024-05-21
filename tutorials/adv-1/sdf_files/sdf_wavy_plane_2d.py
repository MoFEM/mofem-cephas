import numpy as np
from scipy.optimize import brentq

A = 0.0002

wave_len = 2
w = 2. * np.pi / wave_len
ind = 0.00159624

class WavySurface:
    def F(p, x, y, t):
        q = A * (1. - np.cos(w * p)) - ind * t
        lag = 2. * (y - q)
        return 2. * (p - x) - lag * A * w * np.sin(w * p)

    def dist(x, y, a, b, t):
        try:
            p0 = brentq(WavySurface.F, a, b, args=(x, y, t))
            q0 = A * (1. - np.cos(w * p0)) - ind * t
            return np.sqrt((p0 - x)**2 + (q0 - y)**2)
        except ValueError:
            return 1e12

    def sDF(x, y,z, t):
        y0 = A * (1. - np.cos(w * x)) - ind * t
        if np.any(np.abs(y0 - y)) < 1e-9:
            return 0.0

        x = np.abs(x) % wave_len
        x = np.where(x > wave_len / 2, wave_len - x, x)

        eps = 1e-12

        d1 = np.vectorize(WavySurface.dist)(x, y, eps, wave_len / 2 - eps, t) 
        wave_y = A * (1 - np.cos(w * x)) - ind * t
        sgn_d1 = np.sign(wave_y - y)

        wave_y_d2 = -ind * t
        sgn_d2 = np.sign(wave_y_d2 - y)
        d2 = np.sqrt(x**2 + (wave_y_d2 - y)**2) # when xp to be 0
        
        wave_y_d3 = A * (1 - np.cos(w * (wave_len / 2))) - ind * t
        sgn_d3 = np.sign(wave_y_d3 - y)
        d3 = np.sqrt((wave_len / 2 - x)**2 + (A * (1 - np.cos(w * (wave_len / 2))) - ind * t - y)**2) # when xp to be wave_len/2
        
        d_min = np.minimum(d2, d3)
        d_final = np.where(d1 < d_min, d1, np.where(d2 < d3, d2, d3))
        sgn_final = np.where(d1 < d_min, sgn_d1, np.where(d2 < d3, sgn_d2, sgn_d3))
        sdf = sgn_final * d_final
        return sdf


    def gradSDF(x, y, z, t):
        delta = 1e-6
        df_dx = ((WavySurface.sDF(x + delta, y, z, t) - WavySurface.sDF(x - delta, y, z, t)) / (2 * delta)).reshape((-1, 1))
        df_dy = ((WavySurface.sDF(x, y + delta,z, t) - WavySurface.sDF(x, y - delta,z, t)) / (2 * delta)).reshape((-1, 1))
        df_dz = np.zeros_like(df_dy).reshape((-1, 1))
        return np.hstack([df_dx, df_dy, df_dz])

    def hessSDF(x, y, z, t):
        delta = 1e-6
        df2_dx2 = (WavySurface.sDF(x + delta, y,z,t) - 2 * WavySurface.sDF(x, y,z,t) + WavySurface.sDF(x - delta, y,z, t)) / (delta ** 2)
        df2_dy2 = (WavySurface.sDF(x, y + delta,z,t) - 2 * WavySurface.sDF(x, y,z,t) + WavySurface.sDF(x, y - delta,z, t)) / (delta ** 2)
        df2_dxdy = (WavySurface.sDF(x + delta, y + delta,z,t) - WavySurface.sDF(x + delta, y - delta, z, t) - WavySurface.sDF(x - delta, y + delta,z, t) + WavySurface.sDF(x - delta, y - delta,z,t)) / (4 * delta ** 2)

        return np.hstack([df2_dx2.reshape((-1,1)), df2_dxdy.reshape((-1,1)),   np.zeros_like(df2_dy2).reshape((-1,1)),df2_dy2.reshape((-1,1)), np.zeros_like(df2_dy2).reshape((-1,1)), np.zeros_like(df2_dy2).reshape((-1,1))])
    
def sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return WavySurface.sDF(x, y, z, t)

def grad_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return WavySurface.gradSDF(x, y,z, t)

def hess_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return WavySurface.hessSDF(x, y,z, t)