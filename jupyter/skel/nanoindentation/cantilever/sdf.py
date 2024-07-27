import math
import numpy as np

R = 27.5  # radius of the indenter
xc = 0
yc = R
zc = 0
d= 3


# +
def sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return Sphere.sDF(R, xc, yc - d, zc, x, y, z)

def grad_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return Sphere.gradSdf(xc, yc - d, zc, x, y, z)

def hess_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return Sphere.hessSdf(xc, yc - d, zc, x, y, z)


# -

class Sphere:
	def sDF(r, xc, yc, zc, x, y, z):
		return np.sqrt((x - xc)**2 + (y - yc)**2 + (z - zc)**2) - r

	def gradSdf(xc, yc, zc, x, y, z):
		c_val_A = 1./np.sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
		return np.hstack([(c_val_A * (x-xc)).reshape((-1,1)), (c_val_A * (y-yc)).reshape((-1,1)), (c_val_A * (z-zc)).reshape((-1,1))])

	def hessSdf(xc, yc, zc, x, y, z):
		x, y, z = x-xc, y-yc, z-zc
		denom = (x**2 + y**2 + z**2)**(3/2)
		sqrt_denom = np.sqrt(x**2 + y**2 + z**2)
		Hxx = -x**2/denom + 1/sqrt_denom
		Hzx = -x*z/denom
		Hxy = -x*y/denom
		Hyy = -y**2/denom + 1/sqrt_denom
		Hzy = -y*z/denom
		Hzz = -z**2/denom + 1/sqrt_denom
		# xx, yx, zx, yy, zy, zz
		return np.hstack([Hxx.reshape((-1,1)), Hxy.reshape((-1,1)), Hzx.reshape((-1,1)), Hyy.reshape((-1,1)), Hzy.reshape((-1,1)), Hzz.reshape((-1,1))])
