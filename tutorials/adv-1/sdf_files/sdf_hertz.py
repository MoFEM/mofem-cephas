import math
import numpy as np

R = 10 # radius of the indenter
d = 1 # indentation depth

xc = 0
yc = R
zc = 0

def sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return Sphere.sDF(R, xc, yc - d * t, zc, x, y, z)

def grad_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return Sphere.gradSdf(xc, yc - d * t, zc, x, y, z)

def hess_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
    return Sphere.hessSdf(xc, yc - d * t, zc, x, y, z)

class Sphere:
	def sDF(r, xc, yc, zc, x, y, z):
		dx = np.subtract(x, xc)
		dy = np.subtract(y, yc)
		dz = np.subtract(z, zc)
		a = dx**2 + dy**2 + dz**2
		gap = np.sqrt(a) - r
		return gap

	def gradSdf(xc, yc, zc, x, y, z):
		a = (x-xc)**2 + (y-yc)**2 + (z-zc)**2
		c_val = np.sqrt(a)
		c_val_A = 1./c_val
		c_val_dx = c_val_A * (x-xc)
		c_val_dy = c_val_A * (y-yc)
		c_val_dz = c_val_A * (z-zc)

		c_val_dx = c_val_dx.reshape((-1,1))
		c_val_dy = c_val_dy.reshape((-1,1))
		c_val_dz = c_val_dz.reshape((-1,1))
		grad_array = np.hstack([c_val_dx,c_val_dy,c_val_dz])
		return grad_array

	def hessSdf(xc, yc, zc, x, y, z):
		x = x-xc
		y = y-yc
		z = z-zc
		Hxx = -x**2/(x**2 + y**2 + z**2)**(3/2) + 1/np.sqrt(x**2 + y**2 + z**2)
		Hzx = -x*z/(x**2 + y**2 + z**2)**(3/2)
		Hxy = -x*y/(x**2 + y**2 + z**2)**(3/2)
		Hyy = -y**2/(x**2 + y**2 + z**2)**(3/2) + 1/np.sqrt(x**2 + y**2 + z**2)
		Hzy = -y*z/(x**2 + y**2 + z**2)**(3/2)
		Hzz = -z**2/(x**2 + y**2 + z**2)**(3/2) + 1/np.sqrt(x**2 + y**2 + z**2)
		# xx, yx, zx, yy, zy, zz
		Hxx = Hxx.reshape((-1,1))
		Hzx = Hzx.reshape((-1,1))
		Hxy = Hxy.reshape((-1,1))
		Hyy = Hyy.reshape((-1,1))
		Hzy = Hzy.reshape((-1,1))
		Hzz = Hzz.reshape((-1,1))
		hess_array = np.hstack([Hxx, Hxy, Hzx, Hyy, Hzy, Hzz])

		return hess_array
