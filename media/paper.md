---
title: 'MoFEM: An open source, parallel finite element library'
tags:
  - Finite Element
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
  - name: Karol Lewandowski
  - name: Xuan Meng  
  - name: Xiao-Yi Zhou 
  - name: Ignatios Athanasiadis  
  - name: Hoang Nguyen  
  - name: Christophe-Alexandre Chalons-Mouriesse 
  - name: Euan Miur
  - name: Andrei G. Shvarts
  - name: Chris Pearce  
affiliations:
 - name: School of Engineering, University of Glasgow, Glasgow, G12 8QQ
   index: 1
date: 15 April 2019
bibliography: paper.bib
---

# Summary

``MoFEM`` (Mesh Oriented Finite Element Method) is a C++ library supporting the solution of finite elements problems. It is developed to provide free and openß finite element code for engineers, students and academics. ``MoFEM`` is a finite element analysis code tailored for the solution of multi-physics problems with arbitrary levels of approximation, different levels of mesh refinement and optimised for high-performance computing. It is designed to be able to manage complexities related to a heterogeneous order of approximations for L2,H1,H-div and H-curl spaces. ``MoFEM`` is the blend of the Boost MultiIndex containers, MOAB (Mesh Oriented Database)  and PETSc (Portable, Extensible Toolkit for Scientific Computation). ``MoFEM`` is developed in C++ and it is open-source software under the GNU Lesser General Public License.
``MoFEM`` can read and write a number of mesh file formats using functionality provided by MoAB. The current version of MoFEM has full support for CUBIT/TRELIS, TetGEN, Salome/Code_Aster (MED format) for pre-processing and ParaView for post-processing.

# Example users modules

Users modules implemented with ``MoFEM``.
![Example modules.](mofem_modules_examples.png)
Example results of (starting from top left): soap film, solid shell, topology optimisation, phase field fracture, navier stokes, cell traction microscopy, bone remodelling,  configurational fracture, plasticity, mortar contact, magnetostatics, acoustic wave.

# Acknowledgements
![Acknowledgments.](../mofem/doc/figures/Acknowledgments.png)


# References