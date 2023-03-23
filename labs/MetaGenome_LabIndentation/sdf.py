import math

r = 20

#xc = (-58.9 + 12.5)/2
#zc = (-70.5 + 0.9)/2
#yc = -58.9 - r
xc=0
yc=-2.250e+01-r
zc=0.


def sdf(t, x, y, z, tx, ty, tz):
	return Sphere.sDF(r, xc, yc, zc, x, y, z)
	#return CylinderZ.sDF(r, xc, yc, x, y)


def grad_sdf(t, x, y, z, tx, ty, tz):
	return Sphere.gradSdf(xc, yc, zc, x, y, z)
	#return CylinderZ.gradSdf(xc, yc, x, y)


def hess_sdf(t, x, y, z, tx, ty, tz):
	return Sphere.hessSdf(xc, yc, zc, x, y, z)
	#return CylinderZ.hessSdf(xc, yc, x, y)

# Indenters

class yPlane:
	def sDF(shift, xc, yc, x, y):
		return y-shiftl

	def gradSdf(xc, yc, x, y):
		# x, y, z
		return [0, 1, 0]
	
	def hessSdf(xc, yc, x, y):
		# xx, yx, zx, yy, zy, zz
		return [0, 0, 0, 0, 0]


class CylinderZ:
	
	def sDF(r, xc, yc, x, y):
		a = pow(x-xc, 2)+pow(y-yc, 2)
		c_val = math.sqrt(a)-r
		return c_val
	
	def gradSdf(xc, yc, x, y):
		a = pow(x-xc, 2)+pow(y-yc, 2)
		c_val = math.sqrt(a)
		c_val_A = 1/c_val
		c_val_dx = c_val_A * (x-xc)
		c_val_dy = c_val_A * (y-yc)
		# x, y, z
		return [c_val_dx, c_val_dy, 0]
	
	def hessSdf(xc, yc, x, y):
		a = pow(x-xc, 2)+pow(y-yc, 2)
		c_val = math.sqrt(a)
		c_val_A = 1./c_val
		c_val_B = 1./pow(a, 3./2.)
		c_val_dx_dx = c_val_A - c_val_B * pow(x-xc, 2)
		c_val_dx_dy = -c_val_B * (x-xc)*(y-yc)
		c_val_dy_dy = c_val_A - c_val_B * pow(y-yc, 2)
		# xx, yx, zx, yy, zy, zz
		return [c_val_dx_dx, c_val_dx_dy, 0, c_val_dy_dy, 0, 0]

class Sphere:
	
	def sDF(r, xc, yc, zc, x, y, z):
		a = pow(x-xc, 2)+pow(y-yc, 2)+pow(z-zc, 2)
		c_val = math.sqrt(a)-r
		return c_val
	
	def gradSdf(xc, yc, zc, x, y, z):
		a = pow(x-xc, 2)+pow(y-yc, 2)+pow(z-zc, 2)
		c_val = math.sqrt(a)
		c_val_A = 1/c_val
		c_val_dx = c_val_A * (x-xc)
		c_val_dy = c_val_A * (y-yc)
		c_val_dz = c_val_A * (z-zc)
		# x, y, z
		return [c_val_dx, c_val_dy, c_val_dz]
	
	def hessSdf(xc, yc, zc, x, y, z):
		a = pow(x-xc, 2)+pow(y-yc, 2)+pow(z-zc, 2)
		c_val = math.sqrt(a)
		c_val_A = 1./c_val
		c_val_B = 1./pow(a, 3./2.)
		c_val_dx_dx = c_val_A - c_val_B * pow(x-xc, 2)
		c_val_dx_dy = -c_val_B * (x-xc)*(y-yc)
		c_val_dy_dy = c_val_A - c_val_B * pow(y-yc, 2)
		# xx, yx, zx, yy, zy, zz
		return [c_val_dx_dx, c_val_dx_dy, 0, c_val_dy_dy, 0, 0]


