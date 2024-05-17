import numpy as np
from scipy.optimize import brentq

A = 0.2
wave_len = 2
w = 2. * np.pi / wave_len
delta = 1e-12
ind = 0

class WavySurface:
    def F(p, x, y, t):
        q = A * (1 - np.cos(w * p)) - ind * t
        lag = 2. * (y - q)
        return 2. * (p - x) - lag * A * w * np.sin(w * p)

    def dist(x, y, a, b, t):
        try:
            p0 = brentq(WavySurface.F, a, b, args=(x, y, t))
            q0 = A * (1 - np.cos(w * p0)) - ind * t
            return np.sqrt((p0 - x)**2 + (q0 - y)**2)
        except ValueError:
            return 1e12

    def sDF(x, y, t):
        x = np.abs(x) % wave_len
        x = np.where(x > wave_len / 2, wave_len - x, x)

        d1 = WavySurface.dist(x, y, delta, wave_len / 2 - delta, t) 
        wave_y = A * (1 - np.cos(w * x)) - ind * t
        sgn_d1 = np.sign(wave_y - y)

        wave_y_d2 = -ind * t
        sgn_d2 = np.sign(wave_y_d2 - y)
        d2 = np.sqrt((0.0 - x)**2 + (A * (1 - np.cos(w * 0.0)) - ind * t- y)**2)
        
        wave_y_d3 = A * (1 - np.cos(w * (wave_len / 2))) - ind * t
        sgn_d3 = np.sign(wave_y_d3 - y)
        d3 = np.sqrt((wave_len / 2 - x)**2 + (A * (1 - np.cos(w * (wave_len / 2))) - ind * t - y)**2) 
        
        d_min = np.minimum(d2, d3)
        d_final = np.where(d1 < d_min, d1, np.where(d2 < d3, d2, d3))
        sgn_final = np.where(d1 < d_min, sgn_d1, np.where(d2 < d3, sgn_d2, sgn_d3))
        return sgn_final * d_final

    def gradSDF(x, y, t):
        delta = 1e-12
        df_dx = (WavySurface.sDF(x + delta, y, t) - WavySurface.sDF(x - delta, y, t)) / (2 * delta)
        df_dy = (WavySurface.sDF(x, y + delta, t) - WavySurface.sDF(x, y - delta, t)) / (2 * delta)
        return np.array([df_dx, df_dy])

    def hessSDF(x, y, t):
        delta = 1e-6
        df2_dx2 = (WavySurface.sDF(x + delta, y,t) - 2 * WavySurface.sDF(x, y,t) + WavySurface.sDF(x - delta, y, t)) / (delta ** 2)
        df2_dy2 = (WavySurface.sDF(x, y + delta,t) - 2 * WavySurface.sDF(x, y,t) + WavySurface.sDF(x, y - delta, t)) / (delta ** 2)
        df2_dxdy = (WavySurface.sDF(x + delta, y + delta,t) - WavySurface.sDF(x + delta, y - delta, t) - WavySurface.sDF(x - delta, y + delta, t) + WavySurface.sDF(x - delta, y - delta,t)) / (4 * delta ** 2)

        return np.array([[df2_dx2, df2_dxdy], [df2_dxdy, df2_dy2]])
    
def sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return WavySurface.sDF(x, y, t)

def grad_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return WavySurface.gradSDF(x, y, t)

def hess_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return WavySurface.hessSDF(x, y, t)