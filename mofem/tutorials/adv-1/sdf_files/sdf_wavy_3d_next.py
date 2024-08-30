import numpy as np
from scipy.optimize import least_squares

A = 0.01
wave_len = 2.0
w = 2. * np.pi / wave_len
ind = 0.15

class wavy_3D:
    @staticmethod
    def wavy_surface(x, y, t):
        return A * (1.0 - np.cos(w * x) * np.cos(w * y)) - ind * t

    @staticmethod
    def equations(var, x, y, z, t):
        p, q = var
        r = wavy_3D.wavy_surface(p, q, t)
        lag = 2.0 * (z - r)
        eq1 = 2.0 * (p - x) - lag * A * w * np.sin(w * p) * np.cos(w * q)
        eq2 = 2.0 * (q - y) - lag * A * w * np.cos(w * p) * np.sin(w * q)
        return [eq1, eq2]

    @staticmethod
    def sDF(x_arr, y_arr, z_arr, t):
        x_arr = np.abs(x_arr) % wave_len
        x_arr = np.where(x_arr > wave_len / 2.0, wave_len - x_arr, x_arr)
        y_arr = np.abs(y_arr) % wave_len
        y_arr = np.where(y_arr > wave_len / 2.0, wave_len - y_arr, y_arr)
        z_test = wavy_3D.wavy_surface(x_arr, y_arr, t)
        if np.all(np.abs(z_test - z_arr) < 1e-8):
            return np.zeros_like(x_arr)

        sdf_arr = np.zeros_like(x_arr)

        def get_sdf(x, y, z, t, initial_guesses):
            min_sdf = float('inf')
            for initial_guess in initial_guesses:
                solution = least_squares(
                    lambda var: wavy_3D.equations(var, x, y, z, t),
                    initial_guess,
                    bounds=([0.0, 0.0], [wave_len / 2.0, wave_len / 2.0]),
                    verbose=0
                )
                p, q = solution.x
                r = wavy_3D.wavy_surface(p, q, t)
                distance = np.sqrt((p - x)**2 + (q - y)**2 + (r - z)**2)
                sign = np.sign(r - z)
                sdfs = distance * sign
                min_sdf = min(min_sdf, sdfs)
            return min_sdf

        initial_guesses = [[0.0, 0.50], [1.0, 0.50], [0.5, 0.0], [0.5, 1.0]]
        for i in range(len(x_arr)):
            x, y, z = x_arr[i], y_arr[i], z_arr[i]
            sdf_arr[i] = get_sdf(x, y, z, t, initial_guesses)

        return sdf_arr
    @staticmethod
    def gradSDF(x, y, z, t):
        delta = 1e-6
        df_dx = ((wavy_3D.sDF(x + delta, y, z, t) - wavy_3D.sDF(x - delta, y, z, t)) / (2. * delta)).reshape((-1, 1))
        df_dy = ((wavy_3D.sDF(x, y + delta,z, t) - wavy_3D.sDF(x, y - delta,z, t)) / (2. * delta)).reshape((-1, 1))
        df_dz = ((wavy_3D.sDF(x, y,z + delta, t) - wavy_3D.sDF(x, y,z -delta, t)) / (2. * delta)).reshape((-1, 1))
        return np.hstack([df_dx, df_dy, df_dz])

    # symetric hess mat pattern: xx, xy, xz, yy, yz, zz
    @staticmethod
    def hessSDF(x, y, z, t):
        delta = 1e-6
        fxx = (wavy_3D.sDF(x + delta, y, z, t) - 2. * wavy_3D.sDF(x, y, z, t) + wavy_3D.sDF(x - delta, y, z, t)) / (delta ** 2)
        fxy = (wavy_3D.sDF(x + delta, y + delta, z, t) - wavy_3D.sDF(x + delta, y - delta, z, t) - wavy_3D.sDF(x - delta, y + delta, z, t) + wavy_3D.sDF(x - delta, y - delta, z, t)) / (4. * delta ** 2)
        fxz = (wavy_3D.sDF(x + delta, y, z + delta, t) - wavy_3D.sDF(x + delta, y, z - delta, t) - wavy_3D.sDF(x - delta, y, z + delta, t) + wavy_3D.sDF(x - delta, y, z - delta, t)) / (4. * delta ** 2)
        fyy = (wavy_3D.sDF(x, y + delta, z, t) - 2. * wavy_3D.sDF(x, y, z, t) + wavy_3D.sDF(x, y - delta, z, t)) / (delta ** 2)
        fyz = (wavy_3D.sDF(x, y + delta, z + delta, t) - wavy_3D.sDF(x, y + delta, z - delta, t) - wavy_3D.sDF(x, y - delta, z + delta, t) + wavy_3D.sDF(x, y - delta, z - delta, t)) / (4. * delta ** 2)
        fzz = (wavy_3D.sDF(x, y, z + delta, t) - 2. * wavy_3D.sDF(x, y, z, t) + wavy_3D.sDF(x, y, z - delta, t)) / (delta ** 2)
        return np.hstack([fxx.reshape((-1,1)), fxy.reshape((-1,1)), fxz.reshape((-1,1)), fyy.reshape((-1,1)),fyz.reshape((-1,1)),fzz.reshape((-1,1))])
        
def sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return wavy_3D.sDF(x, y, z, t)

def grad_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return wavy_3D.gradSDF(x, y, z, t)

def hess_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return wavy_3D.hessSDF(x, y, z, t)
