SCL-8: Radiation boundary conditions {#jup_radiation_pluto}
=======================================================================
![](https://www.nasa.gov/sites/default/files/thumbnails/image/nh_pluto_10.png)

*Source* [https://www.nasa.gov](https://www.nasa.gov)


# Introduction


In this example, we focus on modelling the *partial differential equation* with a *nonlinear boundary conditions*. We choose the problem of thermal equilibrium of the Pluto planet.The planet model is oversimplified. However, it serves the purpose; this example shows how to linearize partial differential equations and apply the Newton method.  

MoFEM is generic code, and other processes (chemical reactions, radioactive decay, etc.) can be included, to fit better reality. For example, an interesting modelling avenue would be to apply multi-scale analysis, with separate scales for space and time for processes on the surface and subsurface, or model extension to model transport phenomena with reaction-diffusion equation, or near surface convection diffusion problem.

# Part I

## Energy balance


### Conservation law

The universal law of nature, that energy has to be conserved, that with the help of [divergence theorem](https://en.wikipedia.org/wiki/Divergence_theorem) and [principle of locality](https://en.wikipedia.org/wiki/Principle_of_locality), can be expressed for every point of the body as follows,   

\begin{equation}
\sum^2_{i=0} \frac{\partial q_i}{\partial x_i} - f(\mathbf{x}) = 0,\quad\forall \mathbf{x} \in \Omega
\label{eq:conservation1} \tag{1}
\end{equation}

where:
- \f$\mathbf{x}\f$ is point coordinates in [Cartesian coordinate system](https://en.wikipedia.org/wiki/Cartesian_coordinate_system) in body \f$\Omega\f$,
- \f$\mathbf{q}\f$ is energy flux \f$[W \cdot m^{-2}]\f$, and \f$q_i\f$ is coefficient of energy flux, where \f$i=0,1,2\f$.  For example, \f$f(\mathbf{x})\f$ is an energy source, which results from the decay of radioactive materials, or chemical processes, or internal friction caused by changing gravity. We will set source term to zero, to keep it simple. 

Equation (\f$\ref{eq:conservation1}\f$) expressed with the indices applies to the Cartesian coordinate system only, to express it in the more abstract and general form we can write in absolute notation, as follows
\f[
\begin{equation}
\nabla \cdot \mathbf{q} - f(x) = 0,\quad\forall \mathbf{x} \in \Omega
\label{eq:conservation2} \tag{2}
\end{equation}
\f]
where \f$\nabla \cdot\f$ is divergence symbol. Note that we later solve problem in cylindrical coordinate system, were equation expressed \f$(\ref{eq:conservation1})\f$ is not valid for such system, however absolute notation \f$(\ref{eq:conservation2})\f$ holds.

### Thermal conduction

The conection between energy flux and temerature gradient is given by [Fourier's law](https://en.wikipedia.org/wiki/Thermal_conduction#Fourier's_law), 
\f[
\begin{equation}
q_i(\mathbf{x}) = - k(\mathbf{x}) \frac{\partial T(\mathbf{x})}{\partial x_i}
\label{eq:furier1} \tag{3}
\end{equation}
\f]
where:
- \f$T(\mathbf{x})\f$ is temperature at point \f$\mathbf{x}\f$,
- \f$\frac{\partial T(\mathbf{x})}{\partial x_i}\f$ is gradient of temperature, where \f$i=0,1,2\f$,
- \f$k(\mathbf{x})\f$ heat conductivity \f$[W\cdot m^{-2}]\f$ at point \f$\mathbf{x}\f$ in the body.

Furrie's law expresses observation that heat transfer rate is proportional to the negative gradient of temperature, i.e. heat energy flows from hot region to cold region. Furrie's law is a particular expression of the second law of thermodynamics. In the following, we assume that heat conductivity is homogenous through the planet, and its value is the average value for the rock. 

The thermal conduction law expressed in absolute notation takes form as follows
\f[
\begin{equation}
\mathbf{q}(\mathbf{x}) = - k \pmb\nabla T
\label{eq:furier2} \tag{4}
\end{equation}
\f]

## Exchange of energy with the external world

Pluto exchange energy and mass with the universe through the surface, for simplicity, we restrict ourself to two cases, i.e. sun irradiance, planet and outer space radiation. Total power of energy exchange through the surface is given by
\f[
\begin{equation}
P = \int_\Gamma dP = \int_\Gamma dP^S + dP^r \quad [W]
\label{eq:boundary} \tag{5}
\end{equation}
\f]
where:
- \f$\Gamma\f$ is surface of the planet,
- \f$dP^S\f$ is power surface infinitesimal element irradiation received from the sun,
- \f$dP^r\f$ is power surface infinitesimal element radiation/


### Pluto’s insolation


Following [Pluto's insolation history](https://doi.org/10.1016/j.icarus.2014.12.028) we assume that energy flux from sin is 
\f[
\begin{equation}
dP^S(\mathbf{x},t) = I * dA \quad [W]
\label{eq:solar_raddiation1} \tag{6}
\end{equation}
\f]
where:
- \f$I\f$ \f$[W\cdot m^{-2}]\f$ is isolation valuees value take from [Pluto’s insolation history](https://doi.org/10.1016/j.icarus.2014.12.028),
- \f$\mathbf{N}^S(t)\f$ is unit solar ray direction,
- \f$d\mathbf{A}\f$ is surface area element.

<!-- ![](http://mofem.eng.gla.ac.uk/mofem/html/pluto_insolation.png) -->

<img src="http://mofem.eng.gla.ac.uk/mofem/html/pluto_insolation.png"  style="width:100%;"/>
*Pluto's insolation

From [Pluto's insolation history](https://doi.org/10.1016/j.icarus.2014.12.028), Pluto's motion is complex in time due to its elliptical and eccentric orbit. To keep analysis simple, we arbitrary pick data for one full orbit (1767), where north and south pole have averaged isolation values \f$260\f$ \f$[W\cdot m^{-2}]\f$, and on the equator \f$260\f$ \f$[W\cdot m^{-2}]\f$. For simplicity, we will use
\f[
\begin{equation}
I = -0.23 - 0.3 \sin{(\theta)}
\label{eq:solar_raddiation2} \tag{7}
\end{equation}
\f]
where \f$\theta\f$ is altitude of point \f$\mathbf{x}\f$ on Pluto surface, see [Horizontal_coordinate_system](https://en.wikipedia.org/wiki/Horizontal_coordinate_system). 


### Body and space radiation

We will assume that surface of Pluto radiates energy in the space; it also receives energy from space, which has some relict temperature. We assume that outer space temperature is 3K degrees, following Wikipedia (https://en.wikipedia.org/wiki/Outer_space). Balance of energy received by radiation be quantified by [Plank's law](https://en.wikipedia.org/wiki/Planck%27s_law) integrated by all possible wave frequencies, see https://en.wikipedia.org/wiki/Thermal_radiation. Finally, you can express our ideas into the equation
\f[
\begin{equation}
dP^r(\mathbf{x}) = \epsilon \sigma (T^4(\mathbf{x})-T_S^4) dA \quad [W]
\label{eq:solar_raddiation4} \tag{8}
\end{equation}
\f]
where:
- \f$dP^r\f$ is net thermal surface element radiation,
- \f$\epsilon\f$ is emissivity,
- \f$\sigma\f$ is [Stefan–Boltzmann constant](https://en.wikipedia.org/wiki/Stefan–Boltzmann_law),
- \f$T(\mathbf{x})\f$ is temperature on Pluto surface at position \f$\mathbf{x}\f$,
- \f$T_S\f$ is temperature of space,
- \f$d\mathbf{A}\f$ is surface area element.
- \f$\mathbf{x}\f$ is position.

For simplicity, we will set \f$\epsilon=1\f$, assuming that radiation if for the black body. That might be true for outer space radiation but is not true for Pluto surface. However, it could be good enough approximation, considering that it is very cold on this planet, and the surface is rough, porous with the potentially large intrinsic surface area.

## Equations in strong form

Putting all equations together, we get boundary value problem in the indicess for Cartesian coordinate system
\f[
\begin{equation}
\left\{
\begin{array}{l}
\frac{\partial}{\partial x_i} \left( k \frac{\partial T}{\partial x_i} \right)\quad \forall \mathbf{x} \in \Omega, \\
q_i dA_i = dP^S + dP^r \quad\forall \mathbf{} \in \Gamma
\end{array} 
\right.
\label{eq:strong1} \tag{9}
\end{equation}
\f]
Note that we drop \f$\sum\f$ symbol in (\f$\ref{eq:conservation1}\f$), following [Einstein notation](https://en.wikipedia.org/wiki/Einstein_notation). Also, following Einstein notation, \f$q_i dA_i = \sum_{i=0}^2 q_i dA_i\f$.

When problem is expressed in that form, we call it that it is strong form. For clarity, we express above equation in coordinate system independent absolute notation 
\f[
\begin{equation}
\left\{
\begin{array}{l}
\nabla \cdot \left( k \pmb\nabla T \right) = 0\quad\forall \mathbf{x} \in \Omega, \\
\mathbf{q}\cdot d\mathbf{A} = dP^S + dP^r \quad\forall \mathbf{x} \in \Gamma.   
\label{eq:strong2} \tag{10}
\end{array} 
\right.
\end{equation}
\f]
This is a nonlinear partial differential equation since radiation is a nonlinear function of unknown distribution of temperature. You can notice as well that time does not explicitly is involved in the equations. This is a consequence of the assumptions above. We are looking at long periods, like time over several Plutos years (Pluto orbit has ~248 years), where average values are our initial interest.


<!-- 
# Part II

```python
# Running code

# -ts_dt step size
# -ts_max_time 1 final time, should be always 1, since that correspond to accurate solar irradiation
# -snes_rtol relative tolearnce
# -snes_atil absolute tolerance
# -snes_linesearch_type line searcher type (see https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESLineSearchType.html)

import os

user_name=!whoami
user_name=user_name[0]
print(user_name)

if user_name == 'root':
    home_dir = '/mofem_install'
else:
    home_dir=os.environ['HOME'] 

print(home_dir)

dt = 1e-1
snes_rtol = 1e-12
snes_atol = 1e-12
log_file = "log.txt"
wd='%s/um_view/tutorials/scl-8' % home_dir

!cd {wd} && ./radiation \
-file_name pluto.h5m \
-ts_dt {dt} \
-ts_max_time 1 \
-snes_rtol {snes_rtol} \
-snes_atol {snes_atol} \
-snes_max_it 20 \
-snes_linesearch_type bt 2>&1 | tee {log_file}
```

```python
# Convert h5m file into VTK file
!cd {wd} && /mofem_install/um_view/bin/mbconvert out_radiation.h5m out_radiation.vtk
!ls {wd}
```

```python
newton_log_file="newton_log_file.txt"
!cd {wd} && grep "SNES Function norm" {log_file} > {newton_log_file} 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
newton_data=pd.read_fwf(wd+"/"+newton_log_file, header=None)

plt.rcParams['figure.figsize'] = [15, 10]
plt.plot(newton_data[7].to_numpy(),'r^-')
plt.title('Neton method convergence')
plt.ylabel('absolute residial')
plt.xlabel('accumulated iterations')
plt.yscale('log')
plt.grid(True)
plt.show()
```

```python
# Plotting solution on mesh and temperature distribution

import os


# Plot solution
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
pv.start_xvfb()

mesh = pv.read(wd+'/'+'out_radiation.vtk')
my_cmap = plt.cm.get_cmap("jet", 12)

# Take a screen shot
mesh.plot(
    screenshot=wd+'/myscreenshot.png',
    show_grid=True,
    show_edges=True,
    cpos="xy",
    scalars="T", 
    smooth_shading=True, 
    cmap=my_cmap,jupyter_backend='pythreejs')

# Print a screen shotOOO
# import matplotlib.pyplot as plt
# import matplotlib.image as mpimg
# img = mpimg.imread(wd+'/myscreenshot.png')
# imgplot = plt.imshow(img)
# plt.show()
``` 
-->

# Part II

## Equations in weak form

Working with partial differential equations (PDEs) is hard for a computer. Computers are excellent with adding, subtracting and multiplying, but not with solving abstract PDEs. Over many years, many mathematical methods has bee invented which can transform PDEs into the algebraic equations, and then series of linear equations, which can be solved efficiently by computer.

One of the methods which change PDEs into the algebraic equation is the [Galerkin method](https://en.wikipedia.org/wiki/Galerkin_method). If Galerie is used with a particular method of dividing the computational domain and construction of the piecewise polynomial base, it leads to the Finite Element Method (FEM).  Here we skip those details related to FEM and focus on formulation problem in a weak form and linearisation.

We take first equation from the strong form, multiply it by some sufficiently smooth functions \f$v\f$, multiply by it energy conservation equation, and integrate over domain, as result we get
\f[
\begin{equation}
\int_\Omega  v \frac{\partial}{\partial x_i} \left( k \frac{\partial T}{\partial x_i} \right) \textrm{d}V = 0.
\label{eq:weak1} \tag{11}
\end{equation}
\f]
We call function \f$v\f$ testing function, since it test if we satisfy equation \f$(\ref{eq:weak1})\f$. Of curse testing equation using one function \f$v\f$ is not enough, so we can take one by one many different \f$v\f$ functions and check if we fulfil equation \f$(\ref{eq:weak1})\f$ for each of them. In fact, we can have an infinite number of test functions \f$v\f$ against which we can test eq. \f$(\ref{eq:weak1})\f$. However,  that is impractical, since testing for an infinite number of functions will take forever. The solution is the discretisation of the problem, i.e. test functions and temperature can be approximately expressed by a finite number of base functions (e.g. shape functions in finite elements), multiplied by coefficients (e.g. degrees of freedom in finite elements). 

In the next step, we can [integrate by parts](https://en.wikipedia.org/wiki/Integration_by_parts) as follows
\f[
\begin{equation}
-
\int_\Omega  \frac{\partial v}{\partial x_i}   \left( k \frac{\partial T}{\partial x_i} \right) \textrm{d}V 
+
\int_\Gamma v q_i dA_i
= 0.
\label{eq:waak2} \tag{12}
\end{equation}
\f]
and replace the second term using the second equation from the strong form in part I, as result w have
\f[
\begin{equation}
-
\int_\Omega  \frac{\partial v}{\partial x_i}   \left( k \frac{\partial T}{\partial x_i} \right) \textrm{d}V 
+
\int_\Gamma v \left( dP^S + dP^r \right)
= 0
\label{eq:weak4} \tag{13}
\end{equation}
\f]
and finally
\f[
\begin{equation}
r(T(\mathbf{x})) = -
\int_\Omega  \frac{\partial v}{\partial x_i} \left( k \frac{\partial T}{\partial x_i} \right) \textrm{d}V 
+
\int_\Gamma v dP^S
+
\int_\Gamma v \left[
\epsilon \sigma (T^4(\mathbf{x})-T^4_S) 
\right]
dA
= 0
\label{eq:weak5} \tag{14}
\end{equation}
\f]
Note that now we have one equation, instead of two-equation, for body domain and boundary, we have one equation expressed in the integral form which has to be fulfilled for set of test functions \f$v\f$. When the problem is expressed in that way is called a weak form, and boundary condition applied for heat flux is called natural boundary condition. 

You can notice that initially second derivative was applied twice to the field of temperature \f$T(\mathbf{x})\f$, first time to calculate gradient, form which we get heat flux, the second time to divergence. By applying integration by parts, now one derivative form field of temperature \f$T(\mathbf{x})\f$ is "moved" to the test function \f$v(\mathbf{x})\f$. That is important, since not every function can have calculated integral when the second derivative is applied, like in (\f$\ref{eq:weak1}\f$). 

Suppose we like to search for a solution among a set of functions. It is much easier, since it is a bigger set, to search among functions that can have integrated only first derivative, such functions can be piecewise polynomials. 

Note, that piecewise polynomial enables us to approximate functions which have gradients jumps. Such functions have no integrable second derivatives, i.e. are not smooth enough. Moreover, piecewise polynomial is much easier to construct and implement efficiently, in particular in two or three dimensions. 

<img src="https://www.universetoday.com/wp-content/uploads/2008/05/Internal_Structure_of_Pluto.jpg"  style="width:100%;"/>
<!-- ![](https://www.universetoday.com/wp-content/uploads/2008/05/Internal_Structure_of_Pluto.jpg) -->

*Source* www.universetoday.com

Suppose the solution to the problem is a function of temperature, which is [analytical](https://en.wikipedia.org/wiki/Analytic_function). Then a solution of strong form equation and weak form equation is the same in the sense of some appropriate norm. 

However weak form solution is more general and is applicable for problems where temperature distribution is non-analytical. For example, if we like to model Pluto that it has several internal layers, rock core, and outer core layers. Then weak form equation would give correct since at boundary between layers we could observe jumps of the temperature gradient, which can not be represented by an analytical function. From the other hand solution of strong form, equations will demand to divide domains into parts for each layer, and supplement equations by joining each part, enforcing continuity of displacement, and fluxes.  Note that in the model showed above layers are not modelled, and the whole planet is modelled as a homogenous rock.

## Linearisation

You can notice that equation \f$(\ref{eq:weak5})\f$ is nonlinear, residual \f$r(T(\mathbf{x}))\f$, which non-linear function of temperature, i.e. you can find terms like \f$T^4\f$ in the equation. Thus we have to linearise \f$r(T(\mathbf{x}))\f$, with the help comes [Taylor series](https://en.wikipedia.org/wiki/Taylor_series) expansion, which will lead to Newton method. 

We express distribution of temperature as previous value at Newton iteration \f$T_k(\mathbf{x})\f$ and temperature correction  \f$\delta T_{k+1}(\mathbf{x})\f$, such that current temperature is
\f[
\begin{equation}
T_{k+1}(\mathbf{x}) = T_k(\mathbf{x})+\delta T_{k+1}(\mathbf{x}).
\label{eq:taylor1} \tag{15}
\end{equation}
\f]
For convenience, we will write just \f$T_k\f$, instead of \f$T_k(\mathbf{x})\f$. 
Next, using Taylor series expansion, we get
\f[
\begin{equation}
r(T_k+\delta T_{k+1})=
r(T_k) 
+ \left.\frac{\partial r}{\partial T}\right|_{T=T_k} \delta T_{k+1} 
+ \frac{1}{2!}\left.\frac{\partial^2 r}{\partial T^2}\right|_{T=T_k} (\delta T_{k+1})^2
+ \frac{1}{3!}\left.\frac{\partial^3 r}{\partial T^3}\right|_{T=T_k} (\delta T_{k+1})^3
+ \dots
= 0
\label{eq:taylor2} \tag{16}
\end{equation}
\f]
This infinite Taylor series, however, if \f$\delta T_{k+1}\f$ is small, then \f$(\delta T_{k+1})^2\f$ is very small, and following powers are even smaller. Using that observation we can drop all terms beyond the linear term and use \f$\delta T_{k+1}\f$ to correct residual, as follow 
\f[
\begin{equation}
r(T_k+\delta T_{k+1}) \approx 
r(T_k) + \left.\frac{\partial r}{\partial T}\right|_{T=T_k} \delta T_{k+1} = 0. 
\label{eq:taylor3} \tag{17}
\end{equation}
\f]
Solving equation (\f$\ref{eq:taylor2}\f$) is hard, it is much easier to solve (\f$\ref{eq:taylor3}\f$). Solution of (\f$\ref{eq:taylor2}\f$) can be obtained by solving (\f$\ref{eq:taylor3}\f$) iteratively. Starting with \f$k=0\f$, then \f$k=1,2,3,\dots\f$, until \f$r(T_k+\delta T_{k+1})\f$ is small as much we like. If \f$r(T_k+\delta T_{k+1})\f$ is zero, or approaching zero, then solution we have solution of (\f$\ref{eq:weak5}\f$).

This method is called the Newton method. Iteration is repeatedly repeating, the same thing in the loop, and computers are good in such a thing. You can notice that (\f$\ref{eq:taylor3}\f$) is an equation, where field \f$\delta T_{k+1}\f$ is unknown which we will try to find. The computer will do iterations; however, at each iteration, we will have to transform (\f$\ref{eq:taylor3}\f$) into a suitable form for the computer. Next subsection is about that. 

We know distribution of temperature at iteration \f$k\f$, i.e. \f$T_k\f$, since it assumed that it was calculated at previous iteration. We know what is form of \f$r(T_k)\f$, that is shown in eq. \f$(\ref{eq:weak5})\f$, to complete that section we need to derive how second term in (\f$\ref{eq:taylor3}\f$). That can be done simply by calculating derivative only nonlinear term \f$T^4\f$, 
\f[
\begin{equation}
\left.\frac{\partial r}{\partial T}\right|_{T=T_k} \delta T_{k+1} 
= -
\int_\Omega  \frac{\partial v}{\partial x_i} \left( k \frac{\partial \delta T_{k+1}}{\partial x_i} \right) \textrm{d}V 
+
\int_\Gamma v \left[
\epsilon \sigma \left\{ 4T^3_k \right\} \delta T_{k+1}
\right]
dA
\label{eq:taylor4} \tag{18}
\end{equation}
\f]
## Approximation

Next step of the procedure of changing nonlinear PDE into a series of linear equations, which computer can solve. First, for some class of problems, the unknown function of temperature in the body can be expressed by sum of some simpler functions as follows
\f[
\begin{equation}
\begin{split}
T_k(\mathbf{x}) 
&\approx
T^h_k(\mathbf{x}) =
\sum_{a=0}^{N-1} \phi_a(\mathbf{x}) \overline{T}^a_k \\
\delta T_k(\mathbf{x}) 
&\approx
\delta T^h_{k+1}(\mathbf{x}) =
\sum_{a=0}^{N-1} \phi_a(\mathbf{x}) \delta\overline{T}^a_{k+1} \\
v(\mathbf{x}) 
&\approx
v^h(\mathbf{x}) =
\sum_{b=0}^{N-1} \phi_b(\mathbf{x}) \overline{v}^b
\end{split}
\label{eq:approx1} \tag{19}
\end{equation}
\f]
where \f$N\f$ is a number of unknowns coefficients, e.g. \f$\overline{T}^a_{k+1}\f$, and in finite elements it is called number of degrees of freedom. In case of classical Lagrange finite elements, \f$N\f$ in this particular problem is equal to number of nodes, and \f$\phi_a\f$ in nodal function at node \f$a\f$. Nodal functions is picewise polynomial. Superscript \f$(\cdot)^h\f$, indicate that this is discretised solution, which depends on mesh which has some characteristic size of the element \f$h\f$. \f$N\f$ is intrinsically linked to size of mesh \f$h\f$, smaller elements in the mesh more nodes, bigger \f$N\f$.

Not going into the details, by increasing \f$N\f$, i.e. a number of base functions \f$\phi_a\f$, we approximation function \f$T^h\f$ converge to function \f$T\f$, in some norm. Norm is used measure distance two functions from each other. It means that solution of our PDE can be approximated with an arbitrary error, i.e. that distance between \f$T\f$ and \f$T^h\f$ is arbitrary small, as small as which we choose to have, depending on the size of the problem, i.e. \f$N\f$. Bigger \f$N\f$ then functions are closer to each other. Of course big \f$N\f$, means bigger problem, and more equations to solve. 

Replacing \f$T_k(\mathbf{x})\f$ and \f$\delta T_{k+1}(\mathbf{x})\f$, by \f$T^h_k(\mathbf{x})\f$ and \f$\delta T^h_{k+1}(\mathbf{x})\f$
\f[
\begin{equation}
\mathbf{r}(T^h_k(\mathbf{x})) = 
\mathbf{r}^h_a = 
-
\int_\Omega  \frac{\partial \phi_b}{\partial x_i} 
\left( k \frac{\partial T^h_k}{\partial x_i} \right) \textrm{d}V 
+
\int_\Gamma \phi_b dP^S
+
\int_\Gamma \phi_b \left[
\epsilon \sigma ((T^h_k)^4-T^4_s) 
\right]
dA
\label{eq:approx2} \tag{20}
\end{equation}
\f]
and
\f[
\begin{equation}
\mathbf{k}(T^h_k(\mathbf{x})) =
\mathbf{k}^h_{ba}
= -
\int_\Omega  \frac{\partial \phi_b}{\partial x_i} \left( k \frac{\partial \phi_a}{\partial x_i} \right) \textrm{d}V 
+
\int_\Gamma \phi_b \left[
\epsilon \sigma \left\{ 4(T^h_k)^3_k \right\} 
\right]
dA
\label{eq:approx3} \tag{21}
\end{equation}
\f]
where \f$\mathbf{r}^h_a\f$ is vector of size \f$N\f$, and \f$\mathbf{k}^h_{ba}\f$ is matrix \f$N\f$ by \f$N\f$. Having above at hand, we get holy grail of computational methods, i.e. system of linear equations, which now computer can solve,
\f[
\begin{equation}
\mathbf{r}(T^h_k + \delta T^h_{k+1}) \approx
\mathbf{r}^h_k + \mathbf{k}_k^h \delta\overline{T}_{k+1} = 0
\quad
\Rightarrow
\quad
\mathbf{k}_k^h \delta\overline{T}_{k+1} = -\mathbf{r}^h_k
\label{eq:approx4} \tag{22}
\end{equation}
\f]
which we can solve for \f$\delta\overline{T}_{k+1}\f$. Note that \f$\mathbf{k}_k^h\f$ and \f$\mathbf{r}^h_k\f$ for step \f$k\f$, and calculated value of unknown vector \f$\delta\overline{T}_{k+1}\f$  is for step \f$k+1\f$. Equation (\f$\ref{eq:approx4}\f$) is series of linear equation, which can solved efficiently.  


# Part III

## Algorithm for the Newton method

Now we will make loops for the Newton method and time stepping. First, we will start with an internal loop of Newton method:

```python
# max_it - maximal number of iterations
# snes_atol - absolut tolerance
# snes_rtol - relative tolerance

import numpy as np
import matplotlib.pyplot as plt

def newton(res, Kt, T, snes_atol, snes_rtol, max_it):
    
    res0 = res_k = res(T)
    print("It 0: Initial residual r(T) = %3.4e" % abs(res_k))

    arr = [1]
    
    for k in range(0,max_it): 
    
        if abs(res_k) < snes_atol:
            print("Convergeded absolute criterion")
            break
        
        if abs(res_k/res0) < snes_rtol:
            print("Convergeded relative criterion")
            break
    
        dT = - (1./Kt(T))*res_k
        T = T+dT
        res_k = res(T)
        arr = np.append(arr, abs(res_k/res0))
        
        print("It %i: T = %f dT = %f" % (k,T,dT))
        print("It %i: Absolute residual r(T) = %3.4e" % (k, abs(res_k)))
        print("It %i: Relative residual r(T) = %3.4e" % (k, abs(res_k/res0)))
        print()
        
    if abs(res_k) > snes_atol and abs(res_k/res0) > snes_rtol:
        print("Newton method not converged")  
        
    return [T,arr]

def res(T):
    flux = 10
    T_space = 3
    return pow(T,4) - pow(T_space,4) - flux

def Kt(T):
    return 4*pow(T,3)    
     
[Tsol,arr] = newton(res, Kt, 3, 1e-14, 1e-14, 10)  

plt.rcParams['figure.figsize'] = [15, 10]
plt.plot(arr,'r^-')
plt.title('Neton method convergence')
plt.ylabel('relative residial')
plt.xlabel('iteration')
plt.yscale('log')
plt.grid(True)
plt.show()

def cal_log(arr, k):
    return np.log(abs(arr[k+1]-arr[k])/abs(arr[k]-arr[k-1]))

def calc_rate(arr, k):
    return cal_log(arr, k)/cal_log(arr, k-1)

print("Estimated order of convergence: %3.4e" % calc_rate(arr, len(arr)-2))
```

This is how Newton method looks; it is evaluating residual \f$\mathbf{r}\f$, matrix \f$\mathbf{k}\f$, at every iteration.  Next, solve a linear system of equations, and finally, check if residual is small. To check convergence,  criteria to stop Newton iterations are set. One criterion checking absolute residual, other relative residual. In some other cases, criterion on the norm of \f$dT\f$ increment can be set, or change of system energy. The matrix  \f$\mathbf{k}\f$ is often called the tangent stiffness matrix (that term comes from structural analysis).

Here we have a relatively simple function. In general, the key difficulty for more complex problems is to calculate correctly tangent matrix \f$\mathbf{k}\f$. That is essential to achieve a quadratic rate of convergence. The Newton method iterations are costly for big problems, and "quadratic" order of convergence allows to obtain solution fast.  

## Substepping

Newton method is conditionally stable; it converges if the initial guess is close enough to the solution. The initial solution has to be in so-called "convergence ball" if it is outside of it, Newton methods diverge, i.e. not yield a meaning solution. In that case, substepping can help. Another technique improving convergence is a line-searcher method, which is used in Pluto example. However, we will skip the description of line-searcher for another occasion. 

If Pluto is not subjected to sun irradiation, then it is easy to guess the solution, since it is in thermal equilibrium with outer space. In that case, the temperature is uniform over the whole planet and is equal to ambient temperature. We can exploit that, by slowly increasing flux from Sun irradiation until we reach desired real value.  

Algorithms can look as follows:

```python
step_arr = []
for t in np.arange(0.0, 1.25, 0.25):
    print("Psudo time %f" % t)
    
    def res(T):
        flux = t*100
        T_space = 3
        return pow(T,4) - pow(T_space,4) - flux
    
    
    [Tsol,arr] = newton(res, Kt, 3, 1e-12, 1e-12, 10)
    step_arr = np.append(step_arr, arr)
    print()
```

<!--
 ```python 
plt.rcParams['figure.figsize'] = [15, 10]
plt.plot(step_arr,'r^-')
plt.title('Neton method convergence')
plt.ylabel('relative residial')
plt.xlabel('acumilated iterations')
plt.yscale('log')
plt.grid(True)
plt.show()
```
-->

<!-- Note that this plot looks very similar to plot obtained for axisymmetric model of Pluto planet directly calculated in MoFEM. -->

<!--
# Final remarks

Look into the C++ code.

```python

```
-->