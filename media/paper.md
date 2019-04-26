---
title: 'MoFEM: An open source, parallel finite element library'
tags:
  - Finite element method
  - Solid mechanics
  - Fluid mechanics
  - Fracture mechanics
  - Biomechanics
  - C++
  - HPC
authors:
  - name: Lukasz Kaczmarczyk
    orcid: 0000-0002-8468-5435
    affiliation: 1 
  - name: Zahur Ullah
    orcid: 0000-0002-1612-9066
    affiliation: 2
  - name: Karol Lewandowski
    affiliation: 1 
  - name: Xuan Meng  
    affiliation: 1 
  - name: Xiao-Yi Zhou 
    affiliation: 3 
  - name: Ignatios Athanasiadis  
    affiliation: 1 
  - name: Hoang Nguyen  
    affiliation: 1 
  - name: Christophe-Alexandre Chalons-Mouriesse 
    affiliation: 1 
  - name: Euan Miur
    affiliation: 1 
  - name: Andrei G. Shvarts
    affiliation: 1 
  - name: Chris Pearce  
    affiliation: 1 
affiliations:
 - name: School of Engineering, University of Glasgow, Glasgow, G12 8QQ
   index: 1
 - name: School of Mechanical & Aerospace Engineering, Queen's University, Belfast, BT7 1NN
   index: 2
 - name: Department of Bridge Engineering, Tongji University, 1239 Siping Road, Shanghai, 200092, China
   index: 3
date: 25 April 2019
bibliography: paper.bib
---


# Introduction

 `MoFEM` (Mesh-Oriented Finite Element Method) is a C++ library for managing
 complexities related to the finite element method which is a popular numerical
 approach for solving partial differential equations (PDEs) arising in various
 physical problems and engineering applications. `MoFEM` is developed to provide
 free and open source finite element codes, incorporating modern approximation
 approaches and data structures, for engineers, students and academics.

  The need for solutions to increasingly complex problems demands the control
  over numerical errors; otherwise, we will be unable to distinguish
  discretisation artefacts from the real physical phenomena. A brute force
  approach based on a pure *h-adaptivity* which relies on the power of parallel
  computing leads to a low polynomial convergence rate. Therefore, it is
  insufficient to have total control over numerical errors. A more sophisticated
  approach was paved by Ivo Babuska et al. [@babuska1992version], who showed
  that if one could increase at the same time the polynomial order and the mesh
  density, i.e. employ *hp-adaptivity*, the exponential convergence is
  achievable. This has been seen as the 'Holy Grail' of the finite element
  method.

  However, raising the order of approximation comes with a cost where the
  algebraic solver time and the matrix assembly time are increased.
  Unfortunately, there is no universal solution to tackle these two difficulties
  simultaneously. To reduce the algebraic solver time, one way is to use
  multi-grid solvers, which can work more efficiently if a hierarchical
  approximation base is available
  [@ainsworth2003hierarchic,@fuentes2015orientation]. This approach is ideal for
  elliptic problems such as solid elasticity. However, for hyperbolic problems,
  e.g. acoustic wave propagation, the efficiency bottleneck could be in the time
  of the matrix assembly. For that case heterogeneous approximation bases, e.g.
  Bernstein-Bézier base [@ainsworth2011bernstein], allowing for a fast numerical
  integration, could be an optimal solution.

  The control of numerical errors is possible if we can estimate the error to
  drive *hp-adaptivity* algorithm. This error estimator needs to be as much
  efficient as possible, and one possible solution is to use mixed finite
  element formulations, where error estimators become a part of the formulation.
  However, the stability of such elements is an issue, which can be addressed by
  the appropriate use of a combination of $H^1$, $H-curl$, $H-div$ and $L^2$
  spaces. Mixed formulations have other advantages including reduced regularity
  of approximation, or the resulting sparse system of equations, that can be
  exploited by problem-tailored solution algorithms.

  `MoFEM` is designed to provide all solutions discussed above for
  *hp-adaptivity*, enabling rapid implementation of the finite element method
  for solving complex multi-domain, multi-scale and multi-physics engineering
  problems. Moreover, it releases users from programming complexities related to
  the bookkeeping of degrees of freedom (DOFs), finite elements, matrix
  assembly, etc.
  

# Design

  Modern finite element software is an 'ecosystem' managing various complexities
  related to mesh and topology, sparse algebra and approximation, integration
  and dense tensor algebra at the integration point level. Fig.
  \ref{fig:ecosystem} shows different components that are employed in `MoFEM`
  ecosystem. `MoFEM` has not developed and will not develop all these
  capabilities from scratch. Instead `MoFEM` integrates advanced scientific
  computing tools for sparse algebra from
  [PETSc](https://www.mcs.anl.gov/petsc/) (Portable, Extensible Toolkit for
  Scientific Computation) [@petsc-web-page], components for handling mesh and
  topology from [MOAB](https://press3.mcs.anl.gov/sigma/moab-library/)
  (Mesh-Oriented Database) [@tautges_moab:2004] and data structures from [Boost
  libraries](https://www.boost.org) [@boost-web-page]. A further illustration of
  how these packages are utilised in `MoFEM` is shown in Fig. \ref{fig:design}.
  Finally, `MoFEM` core library is developed to manage complexities directly
  related to the finite element method. Therefore, each part of this ecosystem
  has its own design objectives and appropriate programming tools from a
  spectrum of solutions can be selected. Resilience of `MoFEM` ecosystem is
  ensured since the underpinning components have sustainable fundings, dynamic
  and established groups of developers and significant user base.

  !['Ecosystem' of `MoFEM` (Adopted from MoFEM webpage [@MoFEMWebPage]).\label{fig:ecosystem}](ecosystem.png){width=80%}

  <!--  MoFEM makes
  PETSc integral part of code by extending PETSc by DMMOFEM interface (several
  other functions work directly on PETSc objects). MoAB from other hand is
  internal data storage.  -->

  <!--  MoFEM focuses attention on complexities related to finite element
  technology and uses abstractions like field entity, DOF (degree of freedom),
  finite element and problem. -->

  <!-- MoFEM software utilises recent advances
  in the finite element technology and modern data structures, enabling the efficient
  solution of complex, multi-domain, multi-scale and multi-physics problems.  
  -->

  Traditional finite element codes are element-centric meaning the type of an
  element defines the approximation space and base. Therefore, they are not able
  to fully exploit the potential of emerging approximation methods. On the
  contrary, the design of data structures for approximation of field variables
  in `MoFEM`  is independent of the specific finite element, e.g. Lagrangian,
  Nedelec, Rivart-Thomas, since finite element is constructed by a set of lower
  dimension entities on which the approximation fields are defined.
  Consequently, different approximation spaces ($H^1$, $H-curl$, $H-div$, $L^2$)
  can be arbitrarily mixed in a finite element to create new capabilities for
  solving complex problems efficiently. It is worth mentioning that the
  approximation space defines the adjacency of DOFs on entities while the number
  of DOFs on entity is independent on approximation base.

   ![Basic design of `MoFEM` (Adopted from MoFEM webpage [@MoFEMWebPage]).\label{fig:design}](basic_design.png){width=80%}
  
 <!--  Moreover, the base on entity is a trace of the base on element,
  and opposite relation works, base on entity is extruded into element. -->

  `MoFEM` data structures enable easy enrichment of approximation fields and
  modification of base functions, for example, in case of resolving singularity
  at the crack front. Applying this technology, it is effortless to construct
  transition elements between domains with different problem formulation and
  physics, e.g. from two-field mixed formulation to single-field formulation, or
  elements with anisotropic approximation order, e.g. with arbitrary high order
  on surface and arbitrary low order through thickness of solid shells). This
  approach also sets the benchmark in terms of how finite element codes are
  implemented, introducing a concept of user-defined data operators acting on
  fields that are associated with entities (vertices, edges, faces and volumes)
  rather on the finite element directly. Such an approach simplifies code
  writing, testing and validation, making the code resilient to bugs.

 

<!-- 
 ```MoFEM``` is a finite element analysis code tailored for the solution of 
 multi-physics problems with arbitrary levels of approximation, different 
 levels of mesh refinement and optimised for high-performance computing. 

 It is designed to be able to manage complexities related to hierarchical basis 
 functions (Legendre, Lobatto or Jacobi polynomials), providing heterogeneous 
 approximation of an arbitary order for L2, H1, H-div and H-curl spaces. 
 ```MoFEM``` incorporates a blend of 
 [Boost Multi-index Containers](https://www.boost.org/doc/libs/1_62_0/libs/multi_index/doc/index.html), 
 [MOAB](https://press3.mcs.anl.gov/sigma/moab-library/) (Mesh Oriented Database) 
 and [PETSc](https://www.mcs.anl.gov/petsc/) (Portable, Extensible Toolkit 
 for Scientific Computation). 
  -->
 
 <!-- `MoFEM` is developed in C++ and it is an 
 open source software under the 
 [GNU Lesser General Public License](https://www.gnu.org/licenses/lgpl.html). 
 `MoFEM` can read and write a number of mesh file formats using functionality
  provided by [MOAB](https://press3.mcs.anl.gov/sigma/moab-library/). Furthemore,
  it has full support for [CUBIT/TRELIS](https://www.csimsoft.com/trelis.jsp),
  [TetGEN](http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1),
  [Salome-Meca](https://www.code-aster.org/spip.php?article303), 
  [Gmsh](http://gmsh.info) for pre-processing and [ParaView](https:www.paraview.org) for post-processing. -->

# Examples

  `MoFEM` core library provides functionality for developing user modules where
  applications for particular finite elements or problems are implemented. User
  module is an independent repository, private or public and independently
  managed by its owner.

  `MoFEM` was initially created with the financial support of the Royal Academy
   of Engineering and EDF Energy to solve the problem of crack propagation in
   the nuclear graphite [@kaczmarczyk2014three],[@kaczmarczyk2017energy]. Over
   time, the domain of applications expanded to include computational
   homogenisation (DURACOMP EPSRC Project EP/K026925/1),
   [@ullah2019unified,@zhou2017stochastic],[@ullah2017multi] bone remodelling
   and fracture (Kelvin Smith Scholarship), modelling of the gel rheology and
   acoustics problems. Moreover, ```MoFEM``` includes an extensive library of
   example applications such as soap film, solid shell, topology optimisation,
   phase field fracture, Navier-Stokes flow, cell traction microscopy, bone
   remodelling, configurational fracture, plasticity, mortar contact,
   magnetostatics and acoustic wave propagation as shown in Fig.
   \ref{fig:examples}.

  ![Examples of user modules implemented using `MoFEM`.\label{fig:examples}](mofem_modules_examples.png){width=100%}

# Acknowledgements

  `MoFEM` development has been supported by EDF Energy Nuclear Generation Ltd.
  (grant no. 4840360333), The Royal Academy of Engineering (grant no.
  RCSRF1516\2\18), DURACOMP EPSRC Project (EP/K026925/1), and Kelvin Smith
  Scholarship programme at University of Glasgow.

# References