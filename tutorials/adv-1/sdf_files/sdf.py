import math
import numpy as np

# SDF Indenter

# Negative level set represents interior of the indenter. 
# This normal points outside of the indenter.



# Functions for MoFEM
def sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
  return list_indenters[0].sDF(x,y,z)


def grad_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
  return list_indenters[0].gradSdf(x,y,z)


def hess_sdf(delta_t, t, x, y, z, tx, ty, tz, block_id):
  return list_indenters[0].hessSdf(x,y,z)

# Example Indenters

class yPlane:
	def __init__(self,Xc,Yc,Zc,shift):
		# Initial Centroid
		self.Xc = Xc
		self.Yc = Yc
		self.Zc = Zc

		# Current Centroid
		self.xc = Xc
		self.yc = Yc
		self.zc = Zc

		# Indenter Dimensions
		self.shift = shift


	def sDF(self, x, y, z):
		return np.subtract(y,self.shift)

	def gradSdf(self, x, y, z):
		dx = np.zeros_like(x).reshape((-1, 1))
		dy = np.ones_like(y).reshape((-1, 1))
		dz = np.zeros_like(z).reshape((-1, 1))
		return np.hstack([dx, dy, dz])
	
	def hessSdf(self, x, y, z):
		zeros = np.zeros_like(x).reshape((-1, 1))
		return np.hstack([zeros for _ in range(6)])  # xx, yx, zx, yy, zy, zz

class CylinderZ:
	def __init__(self, Xc, Yc, Zc, diameter):
		# Initial Centroid
		self.Xc = Xc
		self.Yc = Yc
		self.Zc = Zc

		# Current Centroid
		self.xc = Xc
		self.yc = Yc
		self.zc = Zc

		# Indenter Dimensions
		self.radius = diameter/2

	def sDF(self, x, y, z):
		return np.sqrt((x - self.xc)**2 + (y - self.yc)**2) - self.radius
	
	def gradSdf(self, x, y, z):
		a = (x-self.xc)**2 + (y-self.yc)**2
		c_val_A = 1./np.sqrt(a)
		c_val_dx = c_val_A * (x-self.xc)
		c_val_dy = c_val_A * (y-self.yc)
		c_val_dz = np.zeros_like(c_val_dy)
    	# x, y, z
		return np.hstack([c_val_dx.reshape((-1,1)), c_val_dy.reshape((-1,1)), c_val_dz.reshape((-1,1))])
	
	def hessSdf(self, x, y, z):
		a = (x-self.xc)**2 + (y-self.yc)**2
		c_val_A = 1./np.sqrt(a)
		c_val_B = 1./(a**(3./2.))
		Hxx = c_val_A - c_val_B * (x-self.xc)**2
		Hxy = -c_val_B * (x-self.xc)*(y-self.yc)
		Hyy = c_val_A - c_val_B * (y-self.yc)**2
		zeros = np.zeros_like(Hxx).reshape((-1,1))
    	# Hxx, Hxy, Hzx, Hyy, Hzy, Hzz
		return np.hstack([Hxx.reshape((-1,1)), Hxy.reshape((-1,1)), zeros, Hyy.reshape((-1,1)), zeros, zeros])

class Sphere:
	def __init__(self, Xc, Yc, Zc, diameter):
		# Initial Centroid
		self.Xc = Xc
		self.Yc = Yc
		self.Zc = Zc

		# Current Centroid
		self.xc = Xc
		self.yc = Yc
		self.zc = Zc

		# Indenter Dimensions
		self.radius = diameter/2
	
	def sDF(self, x, y, z):
		return np.sqrt((x - self.xc)**2 + (y - self.yc)**2 + (z - self.zc)**2) - self.radius
	
	def gradSdf(self, x, y, z):
		a = (x-self.xc)**2 + (y-self.yc)**2 + (z-self.zc)**2
		c_val_A = 1./np.sqrt(a)
		c_val_dx = c_val_A * (x-self.xc)
		c_val_dy = c_val_A * (y-self.yc)
		c_val_dz = c_val_A * (z-self.zc)
    	# x, y, z
		return np.hstack([c_val_dx.reshape((-1,1)), c_val_dy.reshape((-1,1)), c_val_dz.reshape((-1,1))])
	
	def hessSdf(self, x, y, z):
		x, y, z = x-self.xc, y-self.yc, z-self.zc
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
  

# Define Indenters below
r = 1

list_indenters = []
list_indenters.append(CylinderZ(0.0,-0.5-r, 0.0, 2*r))
  