 MIX-2: Mixed  formulation for incompressible elasticity {#incompressible_elasticity}
=======================================================================
# Part 1: Theory 

## Mixed  formulation for incompressible elasticity

The governing equations and boundary conditions for the linear elastic problem are:
\f[
\begin{align}
  -\textrm{div} \boldsymbol{\sigma(u)} = \boldsymbol{f} & \quad \textrm{in}\; \Omega \label{eq:momentum} \\
   \boldsymbol{u} = \boldsymbol{\bar{u}} & \quad \textrm{on}\; \partial\Omega_u, \label{eq:dirichlet}\\
   \boldsymbol{\sigma \cdot n} = \boldsymbol{t}  & \quad \textrm{on}\; \partial\Omega_n, \label{eq:neumann}
\end{align}
\f]
With Hooke's constitutive equation defined as:

\begin{align}
\label{eq:constitutive}\boldsymbol{\sigma} = \lambda tr(\boldsymbol{\varepsilon(u)})\boldsymbol{I} + 2 \mu\boldsymbol{\varepsilon(u)}
\end{align}

Where \f$\boldsymbol{\sigma}\f$ is the Cauchy stress, \f$\boldsymbol{\varepsilon}\f$ is the strain, \f$\boldsymbol{\varepsilon(u)}= \frac{1}{2}(\nabla\boldsymbol{u}+[\nabla\boldsymbol{u}]^T)\f$. The Lame constants \f$\lambda\f$ and \f$\mu\f$ are related through Young's modulus \f$E\f$ and the Poisson ratio \f$\nu\f$:

\f[
\begin{align}
\label{eq:lamelambda}\lambda = \frac{\nu E}{(1 + \nu)(1 - 2\nu)}
\end{align}
\f]

\f[
\begin{align}
\label{eq:lamemu}\mu = \frac{E}{2(1 + \nu)}
\end{align}
\f]

The problem above of compressible elasticity is solved using finite elements. However, for the incompressible state (where \f$\nu\f$ = 0.5), a mixed formulation of the linear elastic problem can be considered through the introduction of a new unknown for pressure, \f$p\f$. With \f$p = \lambda \textrm{div}\boldsymbol{u}\f$, the constitutive relation is now defined as:

\f[
\begin{align}
\label{eq:incompressible}\boldsymbol{\sigma} = 2 \mu\boldsymbol{\varepsilon} + \lambda \textrm{div}\boldsymbol{u}\boldsymbol{I} = 2 \mu\boldsymbol{\varepsilon} + p\boldsymbol{I}
\end{align}
\f]

Therefore, we add the following equation to the strong form of the problem:

\f[
\begin{align}
\label{eq:incompressiblestrong}-\textrm{div} \boldsymbol{u} + \frac{1}{\lambda} p = 0 & \quad \textrm{in}\; \Omega
\end{align}
\f]

As standard for Galerkin finite elements, we multiply Eqs. \eqref{eq:momentum} and \eqref{eq:incompressiblestrong} by the appropriate test functions \f$\delta\mathbf{u}\f$ and \f$\delta p\f$ and integrate over the domain \f$\Omega\f$:

\f[
\begin{align}
\label{eq:momentum_int} 2\mu\int_{\Omega}\mathbf{\varepsilon (u)}:\mathbf{\varepsilon(\delta u)} \,\textrm{d}\Omega + \int_{\Omega}p \,\textrm{div}\,\delta\mathbf{u} \,\textrm{d}\Omega &= \int_{\Omega}\delta\mathbf{u}\cdot\mathbf{f}\,\textrm{d}\Omega + \int_{\partial\Omega}\delta\mathbf{u}\cdot\mathbf{t}\,\textrm{d}\partial\Omega \\
\int_{\Omega}\textrm{div}\,\mathbf{u}\,\delta p \,\textrm{d}\Omega -\int_{\Omega}\frac{1}{\lambda} p\, \delta p\,\textrm{d}\Omega  &= 0 
\end{align}
\f]

As a result, we arrive at the following global linear system of equations:

\f[
\begin{align}
\begin{bmatrix}
{\bf{A}} & {\bf{C}} \\
{\bf{C^T}} & {\bf{B}}
\end{bmatrix}
\begin{bmatrix}
\bf{u} \\
\bf{p} \\
\end{bmatrix}
= -
\begin{bmatrix}
{\bf{F}} \\
{\bf{0}} \\
\end{bmatrix}
\end{align}
\f]

where 

\f[
\begin{align}
\mathbf{A} = 2\mu\int_{\Omega}\mathbf{\varepsilon (u)}:\mathbf{\varepsilon(\delta u)}\,\textrm{d}\Omega 
\end{align}
\f]
\f[
\begin{align}
\mathbf{C^T} = \int_{\Omega}\textrm{div}\,\mathbf{ u}\,\delta p\,\textrm{d}\Omega 
\end{align}
\f]
\f[
\begin{align}
\mathbf{B} = -\int_{\Omega}\frac{1}{\lambda} p\, \delta p\,\textrm{d}\Omega 
\end{align}
\f]
\f[
\begin{align}
\mathbf{F} = \int_{\Omega}\delta\mathbf{u}\cdot\mathbf{f}\,\textrm{d}\Omega + \int_{\partial\Omega}\delta\mathbf{u}\cdot\mathbf{t}\,\textrm{d}\partial\Omega
\end{align}
\f]

In the case of incompressible material, the lower right diagonal block of our matrix, \f${\bf{B}} = 0\f$, resulting in the well-known saddle point problem. Careful choice has to be made for the finite element approximations of the displacement and pressure fields to ensure that the discrete problem has a valid and stable solution, i.e. the inf-sup conditions have to be satisfied. 

One element in particular that can be adopted is the well-known Taylor-Hood element, whereby both displacement and pressure fields are continuous and are approximated in Sobolev \f$H^1\f$ functional space. One must note that the polynomial order for the approximation of the displacement field must be greater by one than that of the pressure field to achieve inf-sup stability. Such choice of the element is often denoted as \f$P_k - P_{k-1}\f$ for \f$k \geq 2\f$. See [[Boffi, 2013]](#boffi_2013) for more details on Taylor-Hood elements.


Another element that can be utilised is the Crouzeix-Raviart element. With this element the displacement is approximated in Sobolev 
\f$H^1\f$ functional space and the pressure field is discontinuous in \f$L^2\f$ space. Two order difference is required between displacement and pressure fields for stability, i.e. such element is denoted as \f$P_k - P_{k-2}\f$. For the 2D case, \f$k \geq 2\f$ is sufficient but for 3D where we require degrees of freedom on the faces, \f$k \geq 3\f$.

### References: 

<a id='boffi_2013'></a> 
**[Boffi, 2013]** Boffi D, Brezzi F, Fortin M. Mixed finite element methods and applications, Springer, 2013

<!-- 
# Part 2: Example - Cook's membrane

To demonstrate the discussed above, implementation of mixed elastic problem is applied to the Cook's membrane problem. Both the Taylor-Hood and Crouzeix-Raviart elements are considered as separate analyses for comparison.

In this problem the beam is subjected to a boundary force on the right hand side edge of the body. This is a typical example in solid mechanics for nearly incompressible and compressible analysis. 

<div>
<img src="workshop_cooks_membrane-3.png" width="600">
<a id='fig_1'></a> 
    <center><b>Fig. 1. Cook's membrane.</b></center> -->

