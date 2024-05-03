import math
import numpy as np

d = 0.001596239 # indentation depth
 
xc = 0
yc = 0
zc = 0

def sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return yPlane.sDF(xc, yc - d * t, zc, x, y, z)

def grad_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return yPlane.gradSdf(xc, yc - d * t, zc, x, y, z)

def hess_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return yPlane.hessSdf(xc, yc - d * t, zc, x, y, z)

class yPlane:
	def sDF(xc, yc, zc, x, y, z):
		return np.subtract(yc, y)

	def gradSdf(xc, yc, zc, x, y, z):
		dx = np.zeros_like(x).reshape((-1, 1))
		dy = np.full_like(y, -1).reshape((-1, 1))
		dz = np.zeros_like(z).reshape((-1, 1))
		return np.hstack([dx, dy, dz])
	
	def hessSdf(xc, yc, zc, x, y, z):
		zeros = np.zeros_like(x).reshape((-1, 1))
		return np.hstack([zeros for _ in range(6)])  # xx, yx, zx, yy, zy, zz