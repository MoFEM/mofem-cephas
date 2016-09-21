/** \file UltraWeakTransportElement.hpp
 * \brief Ultra weak implementation of transport element
 *
 */

/*
 * This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#ifndef __ULTRAWEAK_TARNSPORT_ELEMENT_HPP
#define __ULTRAWEAK_TARNSPORT_ELEMENT_HPP

/** \brief Ultra weak transport problem
  \ingroup mofem_ultra_weak_transport_elem

  Note to solve this system you need to use direct solver or proper preconditioner
  for saddle problem.

  It is based on \cite arnold2006differential \cite arnold2012mixed
  <https://www.researchgate.net/profile/Richard_Falk/publication/226454406_Differential_Complexes_and_Stability_of_Finite_Element_Methods_I._The_de_Rham_Complex/links/02e7e5214f0426ff77000000.pdf>

  General problem have form,
  \f[
  \mathbf{A} \boldsymbol\sigma + \textrm{grad}[u] = \mathbf{0} \; \textrm{on} \; \Omega \\
  \textrm{div}[\boldsymbol\sigma] = f \; \textrm{on} \; \Omega
  \f]

*/
struct UltraWeakTransportElement {

  MoFEM::Interface &mField;

  /**
   * \brief definition of volume element

   * It is used to calculate volume integrals. On volume element we set-up
   * operators to cal;ulcerate components of matrix and vector.

   */
  struct MyVolumeFE: public MoFEM::VolumeElementForcesAndSourcesCore {
    MyVolumeFE(MoFEM::Interface &m_field): MoFEM::VolumeElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return 2*order+1; };
  };

  MyVolumeFE feVol;   ///> Instance of volume element

  /** \brief definition of surface element

    * It is used to calculate surface integrals. On volume element are operators
    * evaluating natural boundary conditions.

    */
  struct MyTriFE: public MoFEM::FaceElementForcesAndSourcesCore {
    MyTriFE(MoFEM::Interface &m_field): MoFEM::FaceElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return 2*order+1; };
  };

  MyTriFE feTriFluxValue;   ///< Instance of surface element

  /**
   * \brief construction of this data structure
   */
  UltraWeakTransportElement(MoFEM::Interface &m_field):
  mField(m_field),
  feVol(m_field),
  feTriFluxValue(m_field) {};

  /**
   * \brief destructor
   */
  virtual ~UltraWeakTransportElement() {}

  VectorDouble valuesAtGaussPts;                          ///< values at integration points on element
  ublas::vector<VectorDouble > valuesGradientAtGaussPts;  ///< gradients at integration points on element
  VectorDouble divergenceAtGaussPts;                      ///< divergence at integration points on element
  ublas::vector<VectorDouble > fluxesAtGaussPts;          ///< fluxes at integration points on element

  set<int> bcIndices;
  PetscErrorCode getDirichletBCIndices(IS *is) {
    PetscFunctionBegin;
    std::vector<int> ids;
    PetscErrorCode ierr;
    ids.insert(ids.begin(),bcIndices.begin(),bcIndices.end());
    IS is_local;
    ierr = ISCreateGeneral(
      mField.get_comm(),ids.size(),ids.empty()?PETSC_NULL:&ids[0],PETSC_COPY_VALUES,&is_local
    ); CHKERRQ(ierr);
    ierr = ISAllGather(is_local,is); CHKERRQ(ierr);
    ierr = ISDestroy(&is_local); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

   /**
    * \brief set source term
    * @param  ent  handle to entity on which function is evaluated
    * @param  x    coord
    * @param  y    coord
    * @param  z    coord
    * @param  flux reference to source term set by function
    * @return      error code
    */
  virtual PetscErrorCode getFlux(
    const EntityHandle ent,
    const double x,const double y,const double z,
    double &flux) {
    PetscFunctionBegin;
    flux  = 0;
    PetscFunctionReturn(0);
  }

   /**
    * \brief natural (Dirihlet) boundary conditions (set values)
    * @param  ent   handle to finite element entity
    * @param  x     coord
    * @param  y     coord
    * @param  z     coord
    * @param  value reference to value set by function
    * @return       error code
    */
  virtual PetscErrorCode getResistivity(
    const EntityHandle ent,
    const double x,const double y,const double z,
    ublas::matrix<FieldData> &inv_k
  ) {
    PetscFunctionBegin;
    inv_k.clear();
    for(int dd = 0;dd<3;dd++) {
      inv_k(dd,dd) = 1;
    }
    PetscFunctionReturn(0);
  }

  /**
   * \brief evaluate natural (Dirichlet) boundary conditions
   * @param  ent   entity on which bc is evaluated
   * @param  x     coordinate
   * @param  y     coordinate
   * @param  z     coordinate
   * @param  value vale
   * @return       error code
   */
  virtual PetscErrorCode getBcOnValues(
    const EntityHandle ent,
    const double x,const double y,const double z,
    double &value) {
    PetscFunctionBegin;
    value = 0;
    PetscFunctionReturn(0);
  }

  /**
   * \brief essential (Neumann) boundary condition (set fluxes)
   * @param  ent  handle to finite element entity
   * @param  x    coord
   * @param  y    coord
   * @param  z    coord
   * @param  flux reference to flux which is set by function
   * @return      [description]
   */
  virtual PetscErrorCode getBcOnFluxes(
    const EntityHandle ent,
    const double x,const double y,const double z,
    double &flux) {
    PetscFunctionBegin;
    flux = 0;
    PetscFunctionReturn(0);
  }

  /** \brief data for evaluation of het conductivity and heat capacity elements
    * \infroup mofem_thermal_elem
    */
  struct BlockData {
    double cOnductivity;
    double cApacity;
    Range tEts; ///< constatins elements in block set
  };
  std::map<int,BlockData> setOfBlocks; ///< maps block set id with appropriate BlockData

  PetscErrorCode addFields(const std::string &values,const std::string &fluxes,const int order) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    //Fields
    ierr = mField.add_field(fluxes,HDIV,1); CHKERRQ(ierr);
    ierr = mField.add_field(values,L2,1); CHKERRQ(ierr);

    //meshset consisting all entities in mesh
    EntityHandle root_set = mField.get_moab().get_root_set();
    //add entities to field
    ierr = mField.add_ents_to_field_by_TETs(root_set,fluxes); CHKERRQ(ierr);
    ierr = mField.add_ents_to_field_by_TETs(root_set,values); CHKERRQ(ierr);
    ierr = mField.set_field_order(root_set,MBTET,fluxes,order+1); CHKERRQ(ierr);
    ierr = mField.set_field_order(root_set,MBTRI,fluxes,order+1); CHKERRQ(ierr);
    ierr = mField.set_field_order(root_set,MBTET,values,order); CHKERRQ(ierr);


    PetscFunctionReturn(0);
  }

  /// \brief add finite elements
  PetscErrorCode addFiniteElements(
    const std::string &fluxes_name,
    const std::string &values_name,
    const std::string mesh_nodals_positions = "MESH_NODE_POSITIONS"
  ) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    ErrorCode rval;

    // Set up volume element operators. Operators are used to calculate components
    // of stiffness matrix & right hand side, in essence are used to do volume integrals over
    // tetrahedral in this case.

    // Define element "ULTRAWEAK". Note that this element will work with fluxes_name and
    // values_name. This reflect bilinear form for the problem
    ierr = mField.add_finite_element("ULTRAWEAK",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ULTRAWEAK",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ULTRAWEAK",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ULTRAWEAK",values_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ULTRAWEAK",values_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK",values_name); CHKERRQ(ierr);

    // In some cases you like to use HO geometry to describe shape of the bode, curved edges and faces, for
    // example body is a sphere. HO geometry is approximated by a field,  which can be hierarchical, so shape of
    // the edges could be given by polynomial of arbitrary order.
    //
    // Check if field "mesh_nodals_positions" is defined, and if it is add that field to data of finite
    // element. MoFEM will use that that to calculate Jacobian as result that geometry in nonlinear.
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK",mesh_nodals_positions); CHKERRQ(ierr);
    }
    // Look for all BLOCKSET which are MAT_THERMALSET, takes entities from those BLOCKSETS
    // and add them to "ULTRAWEAK" finite element. In addition get data form that meshset
    // and set cOnductivity which is used to calculate fluxes from gradients of concentration
    // or gradient of temperature, depending how you interpret variables.
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_THERMALSET,it)) {

      // cerr << *it << endl;

      Mat_Thermal temp_data;
      ierr = it->getAttributeDataStructure(temp_data); CHKERRQ(ierr);
      setOfBlocks[it->getMeshsetId()].cOnductivity = temp_data.data.Conductivity;
      setOfBlocks[it->getMeshsetId()].cApacity = temp_data.data.HeatCapacity;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->getMeshsetId()].tEts,true); CHKERRQ_MOAB(rval);
      ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->getMeshsetId()].tEts,"ULTRAWEAK"); CHKERRQ(ierr);

    }

    // Define element to integrate natural boundary conditions, i.e. set values.
    ierr = mField.add_finite_element("ULTRAWEAK_BCVALUE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ULTRAWEAK_BCVALUE",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ULTRAWEAK_BCVALUE",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_BCVALUE",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_BCVALUE",values_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_BCVALUE",mesh_nodals_positions); CHKERRQ(ierr);
    }

    // Define element to apply essential boundary conditions.
    ierr = mField.add_finite_element("ULTRAWEAK_BCFLUX",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ULTRAWEAK_BCFLUX",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ULTRAWEAK_BCFLUX",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_BCFLUX",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_BCFLUX",values_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_BCFLUX",mesh_nodals_positions); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  /**
   * \brief Build problem
   * @param  ref_level mesh refinement on which mesh problem you like to built.
   * @return           error code
   */
  PetscErrorCode buildProblem(BitRefLevel &ref_level) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    //build field
    ierr = mField.build_fields(); CHKERRQ(ierr);
    // get tetrahedrons which has been build previously and now in so called garbage bit level
    Range done_tets;
    ierr = mField.get_entities_by_type_and_ref_level(
      BitRefLevel().set(BITREFLEVEL_SIZE-1),BitRefLevel().set(),MBTET,done_tets
    ); CHKERRQ(ierr);
    // get tetrahedrons which belong to problem bit level
    Range ref_tets;
    ierr = mField.get_entities_by_type_and_ref_level(
      ref_level,BitRefLevel().set(),MBTET,ref_tets
    ); CHKERRQ(ierr);
    ref_tets = subtract(ref_tets,done_tets);
    ierr = mField.build_finite_elements("ULTRAWEAK",&ref_tets,2); CHKERRQ(ierr);
    // get triangles which has been build previously and now in so called garbage bit level
    Range done_faces;
    ierr = mField.get_entities_by_type_and_ref_level(
      BitRefLevel().set(BITREFLEVEL_SIZE-1),BitRefLevel().set(),MBTRI,done_faces
    ); CHKERRQ(ierr);
    // get triangles which belong to problem bit level
    Range ref_faces;
    ierr = mField.get_entities_by_type_and_ref_level(
      ref_level,BitRefLevel().set(),MBTRI,ref_faces
    ); CHKERRQ(ierr);
    ref_faces = subtract(ref_faces,done_faces);
    //build finite elements structures
    ierr = mField.build_finite_elements("ULTRAWEAK_BCFLUX",&ref_faces,2); CHKERRQ(ierr);
    ierr = mField.build_finite_elements("ULTRAWEAK_BCVALUE",&ref_faces,2); CHKERRQ(ierr);
    //Build adjacencies of degrees of freedom and elements
    ierr = mField.build_adjacencies(ref_level); CHKERRQ(ierr);
    //Define problem
    ierr = mField.add_problem("ULTRAWEAK",MF_ZERO); CHKERRQ(ierr);
    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_set_bit("ULTRAWEAK",ref_level); CHKERRQ(ierr);
    // Add element to problem
    ierr = mField.modify_problem_add_finite_element("ULTRAWEAK","ULTRAWEAK"); CHKERRQ(ierr);
    // Boundary conditions
    ierr = mField.modify_problem_add_finite_element("ULTRAWEAK","ULTRAWEAK_BCFLUX"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ULTRAWEAK","ULTRAWEAK_BCVALUE"); CHKERRQ(ierr);
    //build problem
    ierr = mField.build_problems(); CHKERRQ(ierr);
    //mesh partitioning
    //partition
    ierr = mField.partition_problem("ULTRAWEAK"); CHKERRQ(ierr);
    ierr = mField.partition_finite_elements("ULTRAWEAK"); CHKERRQ(ierr);
    //what are ghost nodes, see Petsc Manual
    ierr = mField.partition_ghost_dofs("ULTRAWEAK"); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  struct OpPostProc: MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {
    moab::Interface &postProcMesh;
    std::vector<EntityHandle> &mapGaussPts;
    OpPostProc(
      moab::Interface &post_proc_mesh,
      std::vector<EntityHandle> &map_gauss_pts
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator("VALUES",UserDataOperator::OPCOL),
    postProcMesh(post_proc_mesh),
    mapGaussPts(map_gauss_pts) {
    }
    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data
    ) {
      MoABErrorCode rval;
      PetscFunctionBegin;
      if(type != MBTET) PetscFunctionReturn(0);
      EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
      Tag th_error_flux;
      rval = getVolumeFE()->mField.get_moab().tag_get_handle("ERRORL2_FLUX",th_error_flux); CHKERRQ_MOAB(rval);
      double* error_flux_ptr;
      rval = getVolumeFE()->mField.get_moab().tag_get_by_ptr(
        th_error_flux,&fe_ent,1,(const void**)&error_flux_ptr
      ); CHKERRQ_MOAB(rval);
      {
        double def_val = 0;
        Tag th_error_flux;
        rval = postProcMesh.tag_get_handle(
          "ERRORL2_FLUX",1,MB_TYPE_DOUBLE,th_error_flux,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val
        ); CHKERRQ_MOAB(rval);
        for(vector<EntityHandle>::iterator vit = mapGaussPts.begin();vit!=mapGaussPts.end();vit++) {
          rval = postProcMesh.tag_set_data(th_error_flux,&*vit,1,error_flux_ptr); CHKERRQ_MOAB(rval);
        }
      }
      PetscFunctionReturn(0);
    }
  };

  /**
   * \brief Post process results
   * @return error code
   */
  PetscErrorCode postProc(const string out_file) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    PostProcVolumeOnRefinedMesh post_proc(mField);
    ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("VALUES"); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("FLUXES"); CHKERRQ(ierr);
    post_proc.getOpPtrVector().push_back(new OpPostProc(post_proc.postProcMesh,post_proc.mapGaussPts));
    ierr = mField.loop_finite_elements("ULTRAWEAK","ULTRAWEAK",post_proc);  CHKERRQ(ierr);
    ierr = post_proc.writeFile(out_file.c_str()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  Vec D,D0,F;
  Mat Aij;

  /// \brief create matrices
  PetscErrorCode createMatrices() {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = mField.MatCreateMPIAIJWithArrays("ULTRAWEAK",&Aij); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost("ULTRAWEAK",COL,&D); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost("ULTRAWEAK",COL,&D0); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost("ULTRAWEAK",ROW,&F); CHKERRQ(ierr);
    ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecZeroEntries(D0); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecZeroEntries(D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /// \brief solve problem
  PetscErrorCode solveProblem() {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr = mField.set_global_ghost_vector("ULTRAWEAK",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);


    // Calculate essential boundary conditions

    // clear essental bc indices, it could have dofs from other mesh refinement
    bcIndices.clear();
    // clear operator, just in case if some other operators are left on this element
    feTriFluxValue.getOpPtrVector().clear();
    // set operator to calculate essential boundary conditions
    feTriFluxValue.getOpPtrVector().push_back(new OpEvaluateBcOnFluxes(*this,"FLUXES",D0));
    ierr = mField.loop_finite_elements("ULTRAWEAK","ULTRAWEAK_BCFLUX",feTriFluxValue); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(D0); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(D0); CHKERRQ(ierr);

    // set operators to calculate matrix and right hand side vectors
    feVol.getOpPtrVector().push_back(new OpFluxDivergenceAtGaussPts(*this,"FLUXES"));
    feVol.getOpPtrVector().push_back(new OpValuesAtGaussPts(*this,"VALUES"));
    feVol.getOpPtrVector().push_back(new OpDivTauU_HdivL2(*this,"FLUXES","VALUES",Aij,F));
    feVol.getOpPtrVector().push_back(new OpL2Source(*this,"VALUES",F));
    feVol.getOpPtrVector().push_back(new OpTauDotSigma_HdivHdiv(*this,"FLUXES",Aij,F));
    feVol.getOpPtrVector().push_back(new OpVDotDivSigma_L2Hdiv(*this,"VALUES","FLUXES",Aij,F));
    ierr = mField.loop_finite_elements("ULTRAWEAK","ULTRAWEAK",feVol); CHKERRQ(ierr);

    // calculate right hand side for natural boundary conditions
    feTriFluxValue.getOpPtrVector().clear();
    feTriFluxValue.getOpPtrVector().push_back(new OpRhsBcOnValues(*this,"FLUXES",F));
    ierr = mField.loop_finite_elements("ULTRAWEAK","ULTRAWEAK_BCVALUE",feTriFluxValue); CHKERRQ(ierr);

    // assemble matrices
    ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

    {
      double nrm2_F;
      ierr = VecNorm(F,NORM_2,&nrm2_F); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"nrm2_F = %6.4e\n",nrm2_F);
    }

    // for ksp solver vector is moved into rhs side
    // for snes it is left ond the left
    ierr = VecScale(F,-1); CHKERRQ(ierr);

    //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);
    //MatView(Aij,PETSC_VIEWER_STDOUT_WORLD);
    // std::string wait;
    //std::cin >> wait;
    IS essential_bc_ids;
    ierr = getDirichletBCIndices(&essential_bc_ids); CHKERRQ(ierr);
    ierr = MatZeroRowsIS(Aij,essential_bc_ids,1,D0,F); CHKERRQ(ierr);
    ierr = ISDestroy(&essential_bc_ids); CHKERRQ(ierr);

    // MatView(Aij,PETSC_VIEWER_DRAW_WORLD);
    // std::cin >> wait;

    // Solve
    KSP solver;
    ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
    ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
    ierr = KSPSetUp(solver); CHKERRQ(ierr);
    ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = KSPDestroy(&solver); CHKERRQ(ierr);

    // copy data form vector on mesh
    ierr = mField.set_global_ghost_vector("ULTRAWEAK",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  /// \brief calculate residual
  PetscErrorCode calculateResidual() {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    //calculate residuals
    feVol.getOpPtrVector().clear();
    feVol.getOpPtrVector().push_back(new OpFluxDivergenceAtGaussPts(*this,"FLUXES"));
    feVol.getOpPtrVector().push_back(new OpValuesAtGaussPts(*this,"VALUES"));
    feVol.getOpPtrVector().push_back(new OpValuesGradientAtGaussPts(*this,"VALUES"));
    feVol.getOpPtrVector().push_back(new OpDivTauU_HdivL2(*this,"FLUXES","VALUES",PETSC_NULL,F));
    feVol.getOpPtrVector().push_back(new OpL2Source(*this,"VALUES",F));
    feVol.getOpPtrVector().push_back(new OpTauDotSigma_HdivHdiv(*this,"FLUXES",PETSC_NULL,F));
    feVol.getOpPtrVector().push_back(new OpVDotDivSigma_L2Hdiv(*this,"VALUES","FLUXES",PETSC_NULL,F));
    // feVol.getOpPtrVector().push_back(new OpError_L2Norm(*this,"VALUES"));
    ierr = mField.loop_finite_elements("ULTRAWEAK","ULTRAWEAK",feVol); CHKERRQ(ierr);
    feTriFluxValue.getOpPtrVector().clear();
    feTriFluxValue.getOpPtrVector().push_back(new OpRhsBcOnValues(*this,"FLUXES",F));
    ierr = mField.loop_finite_elements("ULTRAWEAK","ULTRAWEAK_BCVALUE",feTriFluxValue); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    // ierr = VecAXPY(F,-1.,D0); CHKERRQ(ierr);
    // ierr = MatZeroRowsIS(Aij,essential_bc_ids,1,PETSC_NULL,F); CHKERRQ(ierr);
    {
      std::vector<int> ids;
      ids.insert(ids.begin(),bcIndices.begin(),bcIndices.end());
      std::vector<double> vals(ids.size(),0);
      ierr = VecSetValues(F,ids.size(),&*ids.begin(),&*vals.begin(),INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    }
    {
      double nrm2_F;
      ierr = VecNorm(F,NORM_2,&nrm2_F); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"nrm2_F = %6.4e\n",nrm2_F);
      const double eps = 1e-8;
      if(nrm2_F > eps) {
        //SETERRQ(PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"problem with residual");
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode evaluateError() {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    feVol.getOpPtrVector().clear();
    feVol.getOpPtrVector().push_back(new OpFluxDivergenceAtGaussPts(*this,"FLUXES"));
    feVol.getOpPtrVector().push_back(new OpValuesAtGaussPts(*this,"VALUES"));
    feVol.getOpPtrVector().push_back(new OpValuesGradientAtGaussPts(*this,"VALUES"));
    feVol.getOpPtrVector().push_back(new OpError_L2Norm(*this,"VALUES"));
    ierr = mField.loop_finite_elements("ULTRAWEAK","ULTRAWEAK",feVol,0,mField.getCommSize()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /// \brief destroy matrices
  PetscErrorCode destroyMatrices() {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = MatDestroy(&Aij); CHKERRQ(ierr);
    ierr = VecDestroy(&D); CHKERRQ(ierr);
    ierr = VecDestroy(&D0); CHKERRQ(ierr);
    ierr = VecDestroy(&F); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


  /** \brief Assemble \f$ \int_\mathcal{T} \mathbf{A} \boldsymbol\sigma \cdot \boldsymbol\tau \textrm{d}\mathcal{T} \f$
  */
  struct OpTauDotSigma_HdivHdiv: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Mat Aij;
    Vec F;

    OpTauDotSigma_HdivHdiv(
      UltraWeakTransportElement &ctx,
      const std::string field_name,Mat aij,Vec f
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(
      field_name,
      UserDataOperator::OPROW|UserDataOperator::OPROWCOL
    ),
    cTx(ctx),
    Aij(aij),
    F(f) {}

    MatrixDouble NN,transNN,invK;
    VectorDouble Nf;

    /**
     * \brief Assemble matrix
     * @param  row_side local index of row entity on element
     * @param  col_side local index of col entity on element
     * @param  row_type type of row entity, f.e. MBVERTEX, MBEDGE, or MBTET
     * @param  col_type type of col entity, f.e. MBVERTEX, MBEDGE, or MBTET
     * @param  row_data data for row
     * @param  col_data data for col
     * @return          error code
     */
    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscErrorCode ierr;
      PetscFunctionBegin;
      try {
        if(Aij == PETSC_NULL) PetscFunctionReturn(0);
        if(row_data.getFieldData().size()==0) PetscFunctionReturn(0);
        if(col_data.getFieldData().size()==0) PetscFunctionReturn(0);
        EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
        int nb_row = row_data.getFieldData().size();
        int nb_col = col_data.getFieldData().size();
        NN.resize(nb_row,nb_col,false);
        NN.clear();
        FTensor::Index<'i',3> i;
        FTensor::Index<'j',3> j;
        invK.resize(3,3,false);
        // get access to resistivity data by tensor rank 2
        FTensor::Tensor2<double*,3,3> t_inv_k(
          &invK(0,0),&invK(0,1),&invK(0,2),
          &invK(1,0),&invK(1,1),&invK(1,2),
          &invK(2,0),&invK(2,1),&invK(2,2)
        );
        // get base functions
        FTensor::Tensor1<double*,3> t_n_hdiv_row = row_data.getFTensor1HdivN<3>();
        int nb_gauss_pts = row_data.getHdivN().size1();
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          // get integration weight and multiply by element volume
          double w = getGaussPts()(3,gg)*getVolume();
          // in case that HO geometry is defined that below take into account that
          // edges of element are curved
          if(getHoGaussPtsDetJac().size()>0) {
            w *= getHoGaussPtsDetJac()(gg);
          }
          const double x = getCoordsAtGaussPts()(gg,0);
          const double y = getCoordsAtGaussPts()(gg,1);
          const double z = getCoordsAtGaussPts()(gg,2);
          // calculate receptivity (invers of conductivity)
          ierr = cTx.getResistivity(fe_ent,x,y,z,invK); CHKERRQ(ierr);
          for(int kk = 0;kk!=nb_row;kk++) {
            FTensor::Tensor1<const double*,3> t_n_hdiv_col(
              &col_data.getHdivN(gg)(0,HDIV0),
              &col_data.getHdivN(gg)(0,HDIV1),
              &col_data.getHdivN(gg)(0,HDIV2),3
            );
            for(int ll = 0;ll!=nb_col;ll++) {
              NN(kk,ll) += w*t_n_hdiv_row(i)*t_inv_k(i,j)*t_n_hdiv_col(j);
              ++t_n_hdiv_col;
            }
            ++t_n_hdiv_row;
          }
        }
        // matrix is symmetric, assemble other part
        ierr = MatSetValues(
          Aij,
          nb_row,&row_data.getIndices()[0],
          nb_col,&col_data.getIndices()[0],
          &NN(0,0),ADD_VALUES
        ); CHKERRQ(ierr);
        if(row_side != col_side || row_type != col_type) {
          transNN.resize(nb_col,nb_row);
          noalias(transNN) = trans(NN);
          ierr = MatSetValues(
            Aij,
            nb_col,&col_data.getIndices()[0],
            nb_row,&row_data.getIndices()[0],
            &transNN(0,0),ADD_VALUES
          ); CHKERRQ(ierr);
        }


      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

    /**
     * \brief Assemble matrix
     * @param  side local index of row entity on element
     * @param  type type of row entity, f.e. MBVERTEX, MBEDGE, or MBTET
     * @param  data data for row
     * @return          error code
     */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data
    ) {
      PetscErrorCode ierr;
      PetscFunctionBegin;
      try {
        if(data.getFieldData().size()==0) PetscFunctionReturn(0);
        FTensor::Index<'i',3> i;
        FTensor::Index<'j',3> j;
        invK.resize(3,3,false);
        EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
        int nb_row = data.getFieldData().size();
        Nf.resize(nb_row);
        Nf.clear();
        // get access to resistivity data by tensor rank 2
        FTensor::Tensor2<double*,3,3> t_inv_k(
          &invK(0,0),&invK(0,1),&invK(0,2),
          &invK(1,0),&invK(1,1),&invK(1,2),
          &invK(2,0),&invK(2,1),&invK(2,2)
        );
        // get base functions
        FTensor::Tensor1<double*,3> t_n_hdiv = data.getFTensor1HdivN<3>();
        int nb_gauss_pts = data.getHdivN().size1();
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          double w = getGaussPts()(3,gg)*getVolume();
          if(getHoGaussPtsDetJac().size()>0) {
            w *= getHoGaussPtsDetJac()(gg);
          }
          const double x = getCoordsAtGaussPts()(gg,0);
          const double y = getCoordsAtGaussPts()(gg,1);
          const double z = getCoordsAtGaussPts()(gg,2);
          ierr = cTx.getResistivity(fe_ent,x,y,z,invK); CHKERRQ(ierr);
          FTensor::Tensor1<double*,3> t_flux(
            &cTx.fluxesAtGaussPts[gg][0],
            &cTx.fluxesAtGaussPts[gg][1],
            &cTx.fluxesAtGaussPts[gg][2]
          );
          for(int ll = 0;ll!=nb_row;ll++) {
            Nf[ll] += w*t_n_hdiv(i)*t_inv_k(i,j)*t_flux(j);
            ++t_n_hdiv;
          }
        }
        ierr = VecSetValues(
          F,nb_row,&data.getIndices()[0],&Nf[0],ADD_VALUES
        ); CHKERRQ(ierr);
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

  };

  /** \brief Assemble \f$ \int_\mathcal{T} u \textrm{div}[\boldsymbol\tau] \textrm{d}\mathcal{T} \f$
    */
  struct OpDivTauU_HdivL2: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Mat Aij;
    Vec F;

    OpDivTauU_HdivL2(
      UltraWeakTransportElement &ctx,
      const std::string field_name_row,string field_name_col,Mat _Aij,Vec _F
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(
      field_name_row,field_name_col,
      UserDataOperator::OPROW
    ),
    cTx(ctx),
    Aij(_Aij),
    F(_F) {
      //this operator is not symmetric setting this variable makes element
      //operator to loop over element entities (subsimplicies) without
      //assumption that off-diagonal matrices are symmetric.
      sYmm = false;
    }

    VectorDouble divVec,Nf;

    PetscErrorCode doWork(
      int side,EntityType type,
      DataForcesAndSurcesCore::EntData &data
    ) {
      PetscErrorCode ierr;
      PetscFunctionBegin;
      try {
        if(data.getFieldData().size()==0) PetscFunctionReturn(0);
        int nb_row = data.getIndices().size();
        Nf.resize(nb_row);
        Nf.clear();
        bzero(&*Nf.data().begin(),Nf.data().size()*sizeof(FieldData));
        divVec.resize(data.getHdivN().size2()/3,0);
        if(divVec.size()!=data.getIndices().size()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        int nb_gauss_pts = data.getN().size1();
        int gg = 0;
        for(;gg<nb_gauss_pts;gg++) {
          double w = getGaussPts()(3,gg)*getVolume();
          if(getHoGaussPtsDetJac().size()>0) {
            w *= getHoGaussPtsDetJac()(gg);
          }
          ierr = getDivergenceOfHDivBaseFunctions(side,type,data,gg,divVec); CHKERRQ(ierr);
          noalias(Nf) -= w*divVec*cTx.valuesAtGaussPts[gg];
        }
        ierr = VecSetValues(
          F,nb_row,&data.getIndices()[0],
          &Nf[0],ADD_VALUES
        ); CHKERRQ(ierr);
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

  };

  /** \brief \f$ \int_\mathcal{T} \textrm{div}[\boldsymbol\sigma] v \textrm{d}\mathcal{T} \f$
    */
  struct OpVDotDivSigma_L2Hdiv: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Mat Aij;
    Vec F;

    OpVDotDivSigma_L2Hdiv(
      UltraWeakTransportElement &ctx,
      const std::string field_name_row,string field_name_col,Mat _Aij,Vec _F
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(
      field_name_row,field_name_col,
      UserDataOperator::OPROW|UserDataOperator::OPROWCOL
    ),
    cTx(ctx),
    Aij(_Aij),
    F(_F) {

      //this operator is not symmetric setting this variable makes element
      //operator to loop over element entities without
      //assumption that off-diagonal matrices are symmetric.
      sYmm = false;

    }

    ublas::matrix<FieldData> NN,transNN;
    VectorDouble divVec,Nf;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscErrorCode ierr;
      PetscFunctionBegin;
      try {
        if(Aij == PETSC_NULL) PetscFunctionReturn(0);
        if(row_data.getFieldData().size()==0) PetscFunctionReturn(0);
        if(col_data.getFieldData().size()==0) PetscFunctionReturn(0);
        int nb_row = row_data.getFieldData().size();
        int nb_col = col_data.getFieldData().size();
        NN.resize(nb_row,nb_col);
        NN.clear();
        divVec.resize(nb_col,false);
        int nb_gauss_pts = row_data.getHdivN().size1();
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          double w = getGaussPts()(3,gg)*getVolume();
          if(getHoGaussPtsDetJac().size()>0) {
            w *= getHoGaussPtsDetJac()(gg);
          }
          ierr = getDivergenceOfHDivBaseFunctions(
            col_side,col_type,col_data,gg,divVec
          ); CHKERRQ(ierr);
          noalias(NN) -= w*outer_prod(row_data.getN(gg),divVec);
        }
        ierr = MatSetValues(
          Aij,
          nb_row,&row_data.getIndices()[0],
          nb_col,&col_data.getIndices()[0],
          &NN(0,0),ADD_VALUES
        ); CHKERRQ(ierr);
        transNN.resize(nb_col,nb_row);
        ublas::noalias(transNN) = trans(NN);
        ierr = MatSetValues(
          Aij,
          nb_col,&col_data.getIndices()[0],
          nb_row,&row_data.getIndices()[0],
          &transNN(0,0),ADD_VALUES
        ); CHKERRQ(ierr);
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode doWork(
      int side,EntityType type,
      DataForcesAndSurcesCore::EntData &data
    ) {
      PetscErrorCode ierr;
      PetscFunctionBegin;
      try {
        if(data.getIndices().size()==0) PetscFunctionReturn(0);
        if(data.getIndices().size()!=data.getN().size2()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        int nb_row = data.getIndices().size();
        Nf.resize(nb_row);
        Nf.clear();
        int nb_gauss_pts = data.getN().size1();
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          double w = getGaussPts()(3,gg)*getVolume();
          if(getHoGaussPtsDetJac().size()>0) {
            w *= getHoGaussPtsDetJac()(gg);
          }
          noalias(Nf) -= w*data.getN(gg)*cTx.divergenceAtGaussPts[gg];
        }
        ierr = VecSetValues(
          F,nb_row,&data.getIndices()[0],
          &Nf[0],ADD_VALUES
        ); CHKERRQ(ierr);
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

  };

  /** \brief Calculate source therms, i.e. \f$\int_\mathcal{T} f v \textrm{d}\mathcal{T}\f$
  */
  struct OpL2Source: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Vec F;

    OpL2Source(
      UltraWeakTransportElement &ctx,
      const std::string field_name,Vec _F
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    cTx(ctx),
    F(_F) {}

    VectorDouble Nf;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscErrorCode ierr;
        PetscFunctionBegin;
        try {
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
          int nb_row = data.getFieldData().size();
          Nf.resize(nb_row);
          Nf.clear();
          int nb_gauss_pts = data.getHdivN().size1();
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            double w = getGaussPts()(3,gg)*getVolume();
            if(getHoGaussPtsDetJac().size()>0) {
              w *= getHoGaussPtsDetJac()(gg);
            }
            const double x = getCoordsAtGaussPts()(gg,0);
            const double y = getCoordsAtGaussPts()(gg,1);
            const double z = getCoordsAtGaussPts()(gg,2);
            double flux;
            ierr = cTx.getFlux(fe_ent,x,y,z,flux); CHKERRQ(ierr);
            noalias(Nf) += w*data.getN(gg)*flux;
          }
          ierr = VecSetValues(
            F,nb_row,&data.getIndices()[0],
            &Nf[0],ADD_VALUES
          ); CHKERRQ(ierr);
        } catch (const std::exception& ex) {
          std::ostringstream ss;
          ss << "throw in method: " << ex.what() << std::endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }

        PetscFunctionReturn(0);
      }

    };

  /**
   * \brief calculate \f$ \int_\mathcal{S} {\boldsymbol\tau} \cdot \mathbf{n}u \textrm{d}\mathcal{S} \f$

   * This terms comes from differentiation by parts. Note that in this Dirihlet
   * boundary conditions are natural.

   */
  struct OpRhsBcOnValues: public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Vec F;

    OpRhsBcOnValues(
      UltraWeakTransportElement &ctx,const std::string field_name,Vec _F
    ):
    MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    cTx(ctx),
    F(_F) {}

    VectorDouble Nf;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data
    ) {
      PetscErrorCode ierr;
      PetscFunctionBegin;
      try {
        if(data.getFieldData().size()==0) PetscFunctionReturn(0);
        EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
        Nf.resize(data.getIndices().size());
        Nf.clear();
        int nb_gauss_pts = data.getHdivN().size1();
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          const double x = getCoordsAtGaussPts()(gg,0);
          const double y = getCoordsAtGaussPts()(gg,1);
          const double z = getCoordsAtGaussPts()(gg,2);
          double value;
          ierr = cTx.getBcOnValues(fe_ent,x,y,z,value); CHKERRQ(ierr);
          double w = getGaussPts()(2,gg)*0.5;
          if(getNormalsAtGaussPt().size1() == (unsigned int)nb_gauss_pts) {
            noalias(Nf) += w*prod(data.getHdivN(gg),getNormalsAtGaussPt(gg))*value;
          } else {
            noalias(Nf) += w*prod(data.getHdivN(gg),getNormal())*value;
          }
        }
        ierr = VecSetValues(F,data.getIndices().size(),
        &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }
  };

  /**
   * \brief Evaluate boundary conditions on fluxes.
   *
   * Note that Neumann boundary conditions here are essential. So it is opposite
   * what you find in displacement finite element method.
   *

   * Here we have to solve for degrees of freedom on boundary such base functions
   * approximate flux.

   *
   */
  struct OpEvaluateBcOnFluxes: public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {
    UltraWeakTransportElement &cTx;
    Vec X;
    OpEvaluateBcOnFluxes(
      UltraWeakTransportElement &ctx,const std::string field_name,Vec _X
    ):
    MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    cTx(ctx),
    X(_X) {
    }
    virtual ~OpEvaluateBcOnFluxes() {}
    MatrixDouble NN;
    VectorDouble Nf;
    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscErrorCode ierr;
      PetscFunctionBegin;
      try {
        if(data.getFieldData().size()==0) PetscFunctionReturn(0);
        EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
        int nb_dofs = data.getFieldData().size();
        int nb_gauss_pts = data.getHdivN().size1();
        if(3*nb_dofs!=data.getHdivN().size2()) {
          SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong number of dofs");
        }
        NN.resize(nb_dofs,nb_dofs);
        Nf.resize(nb_dofs);
        // get base functions on rows
        FTensor::Index<'i',3> i;

        // Get normal vector. Note that when higher order geometry is set, then
        // face element could be curved, i.e. normal can be different at each integration
        // point.
        double *normal_ptr;
        if(getNormalsAtGaussPt().size1() == (unsigned int)nb_gauss_pts) {
          // HO geometry
          normal_ptr = &getNormalsAtGaussPt(0)[0];
        } else {
          // Linear geometry, i.e. constant normal on face
          normal_ptr = &getNormal()[0];
        }
        // set tensor from pointer
        FTensor::Tensor1<const double*,3> t_normal(normal_ptr,&normal_ptr[1],&normal_ptr[2],3);

        // get base functions
        FTensor::Tensor1<double*,3> t_n_hdiv_row = data.getFTensor1HdivN<3>();

        NN.clear();
        Nf.clear();
        double nrm2 = 0;

        // loop over integration points
        for(int gg = 0;gg<nb_gauss_pts;gg++) {

          // get integration point coordinates
          const double x = getCoordsAtGaussPts()(gg,0);
          const double y = getCoordsAtGaussPts()(gg,1);
          const double z = getCoordsAtGaussPts()(gg,2);

          // get flux on fece for given element handle and coordinates
          double flux;
          ierr = cTx.getBcOnFluxes(fe_ent,x,y,z,flux); CHKERRQ(ierr);
          // get weight for integration rule
          double w = getGaussPts()(2,gg);
          if(gg == 0) {
            nrm2  = sqrt(t_normal(i)*t_normal(i));
          }

          // set tensor of rank 0 to matrix NN elements
          // loop over base functions on rows and columns
          for(int ll = 0;ll!=nb_dofs;ll++) {
            // get column on shape functions
            FTensor::Tensor1<const double*,3> t_n_hdiv_col(
              &data.getHdivN(gg)(0,HDIV0),
              &data.getHdivN(gg)(0,HDIV1),
              &data.getHdivN(gg)(0,HDIV2),3
            );
            for(int kk = 0;kk<=ll;kk++) {
              NN(ll,kk) += w*t_n_hdiv_row(i)*t_n_hdiv_col(i);
              ++t_n_hdiv_col;
            }
            // right hand side
            Nf[ll] += w*t_n_hdiv_row(i)*t_normal(i)*flux/nrm2;
            ++t_n_hdiv_row;
          }

          // If HO geometry increment t_normal to next integration point
          if(getNormalsAtGaussPt().size1() == (unsigned int)nb_gauss_pts) {
            ++t_normal;
            nrm2  = sqrt(t_normal(i)*t_normal(i));
          }

        }

        // get global dofs indices on element
        cTx.bcIndices.insert(data.getIndices().begin(),data.getIndices().end());

        // factor matrix
        cholesky_decompose(NN);
        // solve local problem
        cholesky_solve(NN,Nf,ublas::lower());

        // set solution to vector
        ierr = VecSetValues(X,data.getIndices().size(),&data.getIndices()[0],&Nf[0],INSERT_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

    PetscFunctionReturn(0);
  }

  };

  /**
   * \brief Calculate values at integration points
   */
  struct OpValuesAtGaussPts: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;

    OpValuesAtGaussPts(
      UltraWeakTransportElement &ctx,
      const std::string field_name
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    cTx(ctx) {}

    virtual ~OpValuesAtGaussPts() {}

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      try {

        if(data.getFieldData().size() == 0)  PetscFunctionReturn(0);

        int nb_gauss_pts = data.getN().size1();
        cTx.valuesAtGaussPts.resize(nb_gauss_pts);
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          cTx.valuesAtGaussPts[gg] = inner_prod( trans(data.getN(gg)), data.getFieldData() );
        }

      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /**
   * \brief Calculate gradients of values at integration points
   */
  struct OpValuesGradientAtGaussPts: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;

    OpValuesGradientAtGaussPts(
      UltraWeakTransportElement &ctx,
      const std::string field_name
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    cTx(ctx) {}
    virtual ~OpValuesGradientAtGaussPts() {}

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {
        if(data.getFieldData().size() == 0)  PetscFunctionReturn(0);
        int nb_gauss_pts = data.getDiffN().size1();
        cTx.valuesGradientAtGaussPts.resize(nb_gauss_pts);
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          cTx.valuesGradientAtGaussPts[gg].resize(3,false);
          noalias(cTx.valuesGradientAtGaussPts[gg]) = prod( trans(data.getDiffN(gg)), data.getFieldData() );
        }
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /**
   * \brief calculate flux at integration point
   */
  struct OpFluxDivergenceAtGaussPts: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;

    OpFluxDivergenceAtGaussPts(
      UltraWeakTransportElement &ctx,
      const std::string field_name
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    cTx(ctx) {}

    virtual ~OpFluxDivergenceAtGaussPts() {}

    VectorDouble divVec;
    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      try {
        if(data.getFieldData().size() == 0)  PetscFunctionReturn(0);
        int nb_gauss_pts = data.getDiffN().size1();
        int nb_dofs = data.getFieldData().size();
        cTx.fluxesAtGaussPts.resize(nb_gauss_pts);
        cTx.divergenceAtGaussPts.resize(nb_gauss_pts);
        if(type == MBTRI && side == 0) {
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            cTx.fluxesAtGaussPts[gg].resize(3,false);
            cTx.fluxesAtGaussPts[gg].clear();
          }
          cTx.divergenceAtGaussPts.clear();
        }
        divVec.resize(nb_dofs);
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          ierr = getDivergenceOfHDivBaseFunctions(side,type,data,gg,divVec); CHKERRQ(ierr);
          cTx.divergenceAtGaussPts[gg] += inner_prod(divVec,data.getFieldData());
          noalias(cTx.fluxesAtGaussPts[gg]) += prod(trans(data.getHdivN(gg)),data.getFieldData());
        }
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

  };

  map<double,EntityHandle> errorMap;

  /** \brief calculate error evaluator
    */
  struct OpError_L2Norm: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;

    OpError_L2Norm(
      UltraWeakTransportElement &ctx,
      const std::string field_name):
      MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
      cTx(ctx) {}
    virtual ~OpError_L2Norm() {}

    VectorDouble deltaFlux;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data
    ) {
      PetscErrorCode ierr;
      ErrorCode rval;
      PetscFunctionBegin;
      try {
        if(type != MBTET) PetscFunctionReturn(0);
        int nb_gauss_pts = data.getN().size1();
        EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
        double def_val = 0;
        Tag th_error_flux;
        rval = cTx.mField.get_moab().tag_get_handle(
          "ERRORL2_FLUX",1,MB_TYPE_DOUBLE,th_error_flux,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val
        ); CHKERRQ_MOAB(rval);
        double* error_flux_ptr;
        rval = cTx.mField.get_moab().tag_get_by_ptr(
          th_error_flux,&fe_ent,1,(const void**)&error_flux_ptr
        ); CHKERRQ_MOAB(rval);
        deltaFlux.resize(3,false);
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          double w = getGaussPts()(3,gg)*getVolume();
          if(getHoGaussPtsDetJac().size()>0) {
            w *= getHoGaussPtsDetJac()(gg);
          }
          const double x = getCoordsAtGaussPts()(gg,0);
          const double y = getCoordsAtGaussPts()(gg,1);
          const double z = getCoordsAtGaussPts()(gg,2);
          double flux;
          ierr = cTx.getFlux(fe_ent,x,y,z,flux); CHKERRQ(ierr);
          if(gg == 0) {
            *error_flux_ptr = 0;
          }
          noalias(deltaFlux) = cTx.fluxesAtGaussPts[gg]+cTx.valuesGradientAtGaussPts[gg];
          *error_flux_ptr += w*( inner_prod(deltaFlux,deltaFlux) );
        }
        if(type == MBTET) {
          *error_flux_ptr = sqrt(*error_flux_ptr);
          cTx.errorMap[*error_flux_ptr] = fe_ent;
        }
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

};

#endif //__ULTRAWEAK_TARNSPORT_ELEMENT_HPP

/***************************************************************************//**
 * \defgroup mofem_ultra_weak_transport_elem Ultra weak transport element
 * \ingroup user_modules
 ******************************************************************************/
