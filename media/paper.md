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
 - name: School of Civil Engineering & Geosciences, Newcastle University, Newcastle upon Tyne, NE1 7RU
   index: 3
date: 25 April 2019
bibliography: paper.bib
---


# Introduction and motivation

 `MoFEM` (Mesh-Oriented Finite Element Method) is a C++ library for managing
 complexities related to the finite element method, a popular numerical approach
 for solving partial differential equations (PDEs) arising in various physical
 problems and engineering applications. `MoFEM` is developed to provide free and
 open source finite element codes, incorporating modern approximation approaches and data structures, for engineers, students and academics.

  The need for solutions to increasingly complex problems demands the control
  over numerical errors; otherwise, we will be unable to distinguish
  discretization artefacts from the real physical phenomena. A brute force
  approach based on a pure *h-adaptivity*, relying on the power of parallel
  computing, is leading to a low polynomial convergence rate and, therefore, is
  insufficient to have total control over numerical errors. A more sophisticated
  approach was paved by Ivo Babuska et al. [@babuska1992version], who showed
  that if one can increase at the same time the polynomial order and the mesh
  density, i.e. employ *hp-adaptivity*, the exponential convergence is possible,
  which is seen as the 'Holy Grail' of the finite element method.

  However, raising the order of approximation comes with a cost, increasing the
  algebraic solver time and the matrix assembly time. Unfortunately, there is no
  universal solution to tackle these two difficulties. To reduce the algebraic solver
  time, one way is to use multi-grid solvers, which can work more efficiently if
  a hierarchical approximation base is available, being ideal for elliptic
  problems such as solid elasticity. However, for hyperbolic problems, e.g.
  acoustic wave propagation, the efficiency bottleneck could be in the time of
  the matrix assembly. For that case heterogeneous approximation bases, e.g.
  Bernstein-Bézier base [@ainsworth2011bernstein], allowing for a fast numerical
  integration, could be an optimal solution.

  The control of numerical errors is possible if we can estimate the error to
  drive *hp-adaptivity* algorithm. This error estimator needs to be as much
  efficient as possible, and one possible solution is to use mix finite element
  formulations, where error evaluators become a part of the formulation.
  However, the stability of such elements is an issue, which can be achieved by
  the appropriate use of a combination of H1, H-curl, H-div and L2 spaces.
  Mix-formulations have other advantages, such as reduced regularity of
  approximation, or the resulting sparse system of equations, which can be
  exploited by problem-tailored solution algorithms.

  `MoFEM` is designed to provide all discussed above solutions for
  *hp-adaptivity*, enabling rapid implementation of the finite element method
  for solving complex multi-domain, multi-scale and multi-physics engineering
  freedom (DOFs), finite elements, matrix assembly, etc.
  

# Basic design of MoFEM

  Modern finite element software is an 'ecosystem' managing various complexities
  related to mesh and topology, sparse algebra and approximation, integration
  and dense tensor algebra at the integration point level. `MoFEM` has not
  developed and will not develop all these capabilities from scratch. Instead,
  `MoFEM` integrates advanced scientific computing tools for sparse algebra from [PETSc](https://www.mcs.anl.gov/petsc/)
  (Portable, Extensible Toolkit for Scientific Computation, components for handling mesh and topology from [MOAB](https://press3.mcs.anl.gov/sigma/moab-library/)
  (Mesh-Oriented Database) and data structures from [Boost libraries](https://www.boost.org). Finally, `MoFEM` core
  library is developed to manage complexities directly related to the finite element method. Therefore, each
  part of this ecosystem has its own design objectives and appropriate programming tools from a
  spectrum of solutions can be selected. Resilience of
  `MoFEM` ecosystem is ensured since the underpinning components have
  sustainable funding, dynamic and established group of developers and
  significant user base.

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

  Traditional finite element codes are element-centric (type of an element
  defines the approximation space and base) and therefore cannot exploit the
  potential of emerging approximation methods. On the contrary, in `MoFEM` the
  design of data structures for approximation of field variables is independent
  of the specific finite element (e.g. Lagrangian, Nedelec, Rivart-Thomas),
  since finite element is constructed by a set of lower dimension entities on
  which the approximation fields are defined. Therefore, different approximation
  spaces (H1, H-curl, H-div, L2) can be arbitrarily mixed in a finite element to
  create new capabilities for solving complex problems efficiently (note that
  the approximation space defines adjacency of DOFs on entities, while the
  number of DOFs on entity is independent on approximation base).

  !['Ecosystem' of `MoFEM`.\label{fig:ecosystem}](ecosystem.png){ width=80% }
  
 <!--  Moreover, the base on entity is a trace of the base on element,
  and opposite relation works, base on entity is extruded into element. -->

  `MoFEM` data structures allow to easily enrich approximation fields or modify
  base functions, e.g. to resolve singularity at a crack front. Applying this
  technology, it is effortless to construct transition elements between domains
  with different problem formulation and physics (e.g. from two-field
  mix-formulation to single-field formulation), or elements with anisotropic
  approximation order (e.g. in solid-shells with arbitrary high-order on shell
  surface and arbitrary low-order through thickness). This approach also sets
  the benchmark in terms of how finite element codes are implemented,
  introducing a concept of user-defined data operators acting on fields that are
  associated with entities (vertices, edges, faces and volumes) rather on the
  finite element directly. Such an approach simplifies code writing, testing and
  validation, making the code resilient to bugs.

  ![Basic design of `MoFEM`.\label{fig:design}](basic_design.png){ width=80% }

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
 
 `MoFEM` is developed in C++ and it is an 
 open source software under the 
 [GNU Lesser General Public License](https://www.gnu.org/licenses/lgpl.html). 
 `MoFEM` can read and write a number of mesh file formats using functionality
  provided by [MOAB](https://press3.mcs.anl.gov/sigma/moab-library/). Furthemore,
  it has full support for [CUBIT/TRELIS](https://www.csimsoft.com/trelis.jsp),
  [TetGEN](http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1),
  [Salome-Meca](https://www.code-aster.org/spip.php?article303), 
  [Gmsh](http://gmsh.info) for pre-processing and [ParaView](https:www.paraview.org) for post-processing.

# Examples of user modules

  MoFEM core library provides functionality for developing user modules where applications for particular finite elements or problems are implemented. User module is an independent repository, private or public and independently managed by its owner.

  `MoFEM` is created with the financial support of the Royal Academy of Engineering and EDF Energy to solve a problem of crack propagation in the nuclear graphite [@kaczmarczyk2014three],[@kaczmarczyk2017energy]. Over the time the domain of applications expanded to include computational homogenisation (DURACOMP EPSRC Project EP/K026925/1), [@ullah2019unified],[@zhou2017stochastic],[@ullah2017multi] bone remodelling and fracture (Kelvin Smith Scholarship), modelling of the gels reology and acoustics problems. Moreover, ```MoFEM``` includes an extensive library of example applications such as soap film, solid shell, topology optimisation, phase field fracture, Navier-Stokes flow, cell traction microscopy, bone remodelling, configurational fracture, plasticity, mortar contact, magnetostatics and acoustic wave propagation, see Figure \ref{fig:examples}.

  ![Examples of user modules implemented with `MoFEM`.\label{fig:examples}](mofem_modules_examples.png){ width=100% }

# Acknowledgements

  `MoFEM` development is supported by EDF Energy Nuclear Generation Ltd 
  (grant no. 4840360333), The Royal Academy of Engineering (grant no. 
  RCSRF1516\2\18), DURACOMP EPSRC Project (EP/K026925/1) and Kelvin Smith
  Scholarship programme at University of Glasgow.

# References