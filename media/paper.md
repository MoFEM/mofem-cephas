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
date: 15 April 2019
bibliography: paper.bib
---


# Introduction
 ```MoFEM``` (Mesh-Oriented Finite Element Method) is a C++ library for managing complexities related to the finite element method, a popular approach for solving partial differential equations (PDEs) arising in various physical problems and engeeniring applications. The main purpose behind the development of ```MoFEM``` is to enable rapid implementation of the finite element method for solving complex multi-physics and multi-domain engeneering-related problems utilizing modern hierarchical and heterogeneous approximation approaches.  
 
 # Summary

 ```MoFEM``` is a finite element analysis code tailored for the solution of multi-physics problems with arbitrary levels of approximation, different levels of mesh refinement and optimised for high-performance computing. It is designed to be able to manage complexities related to hierarchical basis functions (Legendre,Lobatto or Jacoby polynomials), providing heterogeneous approximation of an arbitary order for L2, H1, H-div and H-curl spaces. ```MoFEM``` incorporates a blend of [Boost Multi-index Containers](https://www.boost.org/doc/libs/1_62_0/libs/multi_index/doc/index.html), [MOAB](https://press3.mcs.anl.gov/sigma/moab-library/) (Mesh Oriented Database) and [PETSc](https://www.mcs.anl.gov/petsc/) (Portable, Extensible Toolkit for Scientific Computation). ```MoFEM``` is developed in C++ and it is an open source software under the [GNU Lesser General Public License](https://www.gnu.org/licenses/lgpl.html). ```MoFEM``` can read and write a number of mesh file formats using functionality provided by [MOAB](https://press3.mcs.anl.gov/sigma/moab-library/). Furthemore, it has full support for [CUBIT/TRELIS](https://www.csimsoft.com/trelis.jsp), [TetGEN](http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1), [Salome-Meca](https://www.code-aster.org/spip.php?article303), [Gmsh](http://gmsh.info) for pre-processing and [ParaView](https://www.paraview.org) for post-processing.


```MoFEM``` is created with the financial support of the Royal Academy of Engineering and EDF Energy to solve a problem of crack propagation in the nuclear graphite [@kaczmarczyk2014three],[@kaczmarczyk2017energy]. Over the time the domain of applications expanded to include computational homogenisation (DURACOMP EPSRC Project EP/K026925/1), [@ullah2019unified],[@zhou2017stochastic],[@ullah2017multi] bone remodelling and fracture (Kelvin Smith Scholarship), modelling of the gels reology and acoustics problems. Moreover, ```MoFEM``` includes an extensive library of example applications, such as soap film, solid shell, topology optimisation, phase field fracture, Navier-Stokes flow, cell traction microscopy, bone remodelling, configurational fracture , plasticity, mortar contact and magnetostatics.

# Example users modules

Users modules implemented with ``MoFEM``.
![Example modules.](mofem_modules_examples.png)
Example results of (starting from top left): soap film, solid shell, topology optimisation, phase field fracture, navier stokes, cell traction microscopy, bone remodelling,  configurational fracture, plasticity, mortar contact, magnetostatics, acoustic wave.

# Acknowledgements
![Acknowledgments.](../mofem/doc/figures/Acknowledgments.png)


# References