Applications of MoFEM in research and industry (2021-2023) {#applications}
=======================================================================

MoFEM delivers a software development platform which enhances scientific innovation by providing a flexible and adaptable modelling framework, using novel disruptive approaches to long-standing problems in continuum mechanics and tackling conflicting requirements of accuracy and computational efficiency. This is achieved in MoFEM by developing and adopting state-of-the-art FE technologies, for example: \f$H^1\f$-, \f$H(\text{curl})\f$-, \f$H(\text{div}\f$)- and \f$L^2\f$-conforming finite elements equipped with hierarchical, heterogeneous and anisotropic approximation bases; error-driven hp-adaptivity; mesh topology evolution. In addition, MoFEMâ€™s HPC capabilities are supported by its unique data structures that are capable of handling generic multi-field, multi-physics, multi-scale problems and building tailored composite solvers.

Therefore, MoFEM provides users with an effective tool for solving Partial Differential Equations arising in various fields of Engineering and Applied Physics: solid mechanics, fluid mechanics, soft matter physics, heat transfer, electromagnetism, etc. Furthermore, MoFEM features an extendable modular design: while its open-source core library is developed to manage the complexities of FEM, additional user modules are devoted to particular applications. Such a toolkit-like structure allows for independent development of modules with different repositories, owners and licenses.

[TOC]

# Solid mechanics problems {#solid_mechanics}

## Brittle crack propagation under contact loading {#fracture}

\f$\textbf{I. Athanasiadis}^1, \textbf{A. G. Shvarts}^1, \textbf{Å. Kaczmarczyk}^1,\textbf{K. Lewandowski}^1, \textbf{C. Pearce}^1\f$

\f$^1 \textit{Glasgow Computational Engineering Centre, James Watt School of Engineering, University of Glasgow}\f$

<div style="text-align: justify"> We present the first implicit computational framework for simulating crack propagation along contact interfaces and surfaces under load in three-dimensional bodies. We restrict ourselves to brittle fracture and frictionless contact and focus on numerical challenges associated with the coupling of unilateral constraints emerging from the Griffith's criterion and the contact conditions. The formulation is based on the configurational mechanics framework and is solved using the finite element method. The approach utilises a monolithic Arbitrary Lagrangian-Eulerian formulation permitting simultaneous resolution of crack propagation and unilateral contact constraints. Contact is embedded in the model using the well-known mortar contact formulation. Evolving cracks are explicitly modelled as displacement discontinuities within the mesh. Heterogeneous approximation of arbitrary order is used to discretise spatial displacements, enabling \f$hp\f$-adaptive refinement around the crack front and the contact interfaces traversed by the crack. The result is a holistic approach which handles issues associated with thermodynamic consistency, numerical accuracy and robustness of the computational scheme. Several numerical examples are presented to verify the model formulation and implementation; they also highlight how contact pressure and load applied on surfaces traversed by cracks influence their propagation. The robustness of the approach is validated by comparison of our simulations with industrial experiment involving cracks of complex morphologies propagating along contact interfaces between multiple deformable bodies \cite athanasiadis2023computational.
</div>

<img src="jacobs_1.png" alt="Docker - searching for the required Docker image and Tag" width="80%"/>
<a id='figure_1'></a> 
    <center><b>Figure: (a) experimental results for crack propagation in the graphite brick slice caused by contact with the graphite seal ring and the steel loading collar (not shown), (b) schematic side view showing the brick slice, the seal ring and the loading collar, (c) results of the simulation \cite athanasiadis2023computational </b></center>

# Transport problems




# Multiphysics problems