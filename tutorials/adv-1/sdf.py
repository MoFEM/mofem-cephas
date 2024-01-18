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

	def gradSdf(self,x, y,z):
		dx = np.zeros_like(x)
		dy = np.ones_like(y)
		dz = np.zeros_like(z)
		dx = dx.reshape((-1,1))
		dy = dy.reshape((-1,1))
		dz = dz.reshape((-1,1))
		grad_array = np.hstack([dx, dy, dz])
		return grad_array
	
	def hessSdf(self,x, y,z):
		zeros = np.zeros_like(x)
		zeros = zeros.reshape((-1,1))

		hess_array = np.hstack([zeros, zeros, zeros, zeros, zeros, zeros])

 		# xx, yx, zx, yy, zy, zz
		return hess_array

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
		a = (x-self.xc)**2 + (y-self.yc)**2
		gap = np.sqrt(a)-self.radius

		self.normal = self.gradSdf(x,y,z)
		self.dNormal = self.hessSdf(x,y,z)

		return gap
	
	def gradSdf(self, x, y, z):
		a = (x-self.xc)**2 + (y-self.yc)**2
		c_val = np.sqrt(a)
		c_val_A = 1./c_val
		c_val_dx = c_val_A * (x-self.xc)
		c_val_dy = c_val_A * (y-self.yc)
		c_val_dz = np.zeros_like(c_val_dy)
		# x, y, z
		c_val_dx = c_val_dx.reshape((-1,1))
		c_val_dy = c_val_dy.reshape((-1,1))
		c_val_dz = c_val_dz.reshape((-1,1))
		grad_array = np.hstack([c_val_dx,c_val_dy,c_val_dz])
		return grad_array
	
	def hessSdf(self, x, y, z):
		a = (x-self.xc)**2 + (y-self.yc)**2
		c_val = np.sqrt(a)
		c_val_A = 1./c_val
		c_val_B = 1./(a**(3./2.))
		Hxx = c_val_A - c_val_B * (x-self.xc)**2
		Hxy = -c_val_B * (x-self.xc)*(y-self.yc)
		Hyy = c_val_A - c_val_B * (y-self.yc)**2

		Hxx = Hxx.reshape((-1,1))
		Hzx = np.zeros_like(Hxx)
		Hxy = Hxy.reshape((-1,1))
		Hyy = Hyy.reshape((-1,1))
		Hzy = np.zeros_like(Hxx)
		Hzz = np.zeros_like(Hxx)
		hess_array = np.hstack([Hxx, Hxy, Hzx, Hyy, Hzy, Hzz])

		return hess_array

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
		dx = np.subtract(x, self.xc)
		dy = np.subtract(y, self.yc)
		dz = np.subtract(z, self.zc)
		a = (dx)**2+(dy)**2 + (dz)**2
		gap = np.sqrt(a) - self.radius
		
		return gap
	
	def gradSdf(self,x, y, z):
		a = (x-self.xc)**2 + (y-self.yc)**2 + (z-self.zc)**2
		c_val = np.sqrt(a)
		c_val_A = 1./c_val
		c_val_dx = c_val_A * (x-self.xc)
		c_val_dy = c_val_A * (y-self.yc)
		c_val_dz = c_val_A * (z-self.zc)
		# x, y, z
		#size = np.size(x)
		#grad_array = np.empty([size,3])
		c_val_dx = c_val_dx.reshape((-1,1))
		c_val_dy = c_val_dy.reshape((-1,1))
		c_val_dz = c_val_dz.reshape((-1,1))
		grad_array = np.hstack([c_val_dx,c_val_dy,c_val_dz])
		return grad_array
	
	def hessSdf(self,x, y, z):
		x = x-self.xc
		y = y-self.yc
		z = z-self.zc
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
  


# Define Indenters below
r = 1

list_indenters = []
list_indenters.append(CylinderZ(0.0,-0.5-r, 0.0, 2*r))
  