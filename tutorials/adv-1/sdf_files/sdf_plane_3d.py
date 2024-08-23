import math
import numpy as np

d = 0.00217384 # indentation depth
 
xc = 0
yc = 0
zc = 0

def sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return zPlane.sDF(xc, yc, zc - d * t, x, y, z)

def grad_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return zPlane.gradSdf(xc, yc, zc - d * t, x, y, z)

def hess_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return zPlane.hessSdf(xc, yc, zc - d * t, x, y, z)

class zPlane:
	def sDF(xc, yc, zc, x, y, z):
		return np.subtract(zc, z)

	def gradSdf(xc, yc, zc, x, y, z):
		dx = np.zeros_like(x).reshape((-1, 1))
		dy = np.zeros_like(y).reshape((-1, 1))
		dz = np.full_like(z, -1).reshape((-1, 1))
		return np.hstack([dx, dy, dz])
	
	def hessSdf(xc, yc, zc, x, y, z):
		zeros = np.zeros_like(x).reshape((-1, 1))
		return np.hstack([zeros for _ in range(6)])  # xx, yx, zx, yy, zy, zz