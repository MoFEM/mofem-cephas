/** \file MagneticElement.hpp
 * \brief Implementation of magnetic element
 * \ingroup maxwell_element
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

#ifndef __MAGNETICELEMENT_HPP__
#define __MAGNETICELEMENT_HPP__

/**
 * \brief Implementation of magnetostatic problem (basic Implementation)
 * \ingroup maxwell_element
 *
 *  Look for theory and details here:
 *
 *  \cite ivanyshyn2013computation
 *  <www.hpfem.jku.at/publications/szthesis.pdf>
 *
 */
struct MagneticElement {

  MoFEM::Interface &mField;

  /// \brief  definition of volume element
  struct VolumeFE: public MoFEM::VolumeElementForcesAndSourcesCore {
    VolumeFE(MoFEM::Interface &m_field):
    MoFEM::VolumeElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return 2*order+1; };
  };

  // /// \brief  definition of volume element
  // struct VolumeFEReducedIntegration: public MoFEM::VolumeElementForcesAndSourcesCore {
  //   VolumeFEReducedIntegration(MoFEM::Interface &m_field):
  //   MoFEM::VolumeElementForcesAndSourcesCore(m_field) {}
  //   int getRule(int order) { return 2*order+1; };
  // };

  /** \brief define surface element
    *
    */
  struct TriFE: public MoFEM::FaceElementForcesAndSourcesCore {
    TriFE(MoFEM::Interface &m_field): MoFEM::FaceElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return 2*order+1; };
  };

  MagneticElement(MoFEM::Interface &m_field):
  mField(m_field) {
  }
  virtual ~MagneticElement() {}

  /**
   * \brief data structure storing material constants, model parameters, matrices, etc.
   *
   */
  struct BlockData {

    // field
    const string fieldName;
    const string feName;
    const string feNaturalBCName;

    // material parameters
    double mU;  ///< magnetic constant  N / A2
    double ePsilon; ///< regularization paramater

    // Natural boundary conditions
    Range naturalBc;

    // Essential boundary conditions
    Range essentialBc;

    int oRder; ///< approximation order

    // Petsc data
    DM dM;
    Mat A;
    Vec D,F;

    BlockData():
    fieldName("MAGNETIC_POTENTIAL"),
    feName("MAGNETIC"),
    feNaturalBCName("MAGENTIC_NATURAL_BC"),
    mU(1),
    ePsilon(0.1) {
    }
    ~BlockData() {}
  };

  BlockData blockData;

  /**
   * \brief get natural boundary conditions
   * \ingroup maxwell_element
   * @return      error code
   */
  PetscErrorCode getNaturalBc() {
    MoABErrorCode rval;
    PetscFunctionBegin;
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,bit)) {
      if(bit->getName().compare(0,9,"NATURALBC") == 0) {
        Range faces;
        rval = mField.get_moab().get_entities_by_type(bit->meshset,MBTRI,faces,true); CHKERRQ_MOAB(rval);
        rval = mField.get_moab().get_adjacencies(
          faces,1,true,blockData.naturalBc,moab::Interface::UNION
        ); CHKERRQ_MOAB(rval);
        blockData.naturalBc.merge(faces);
      }
    }
    PetscFunctionReturn(0);
  }

  /**
   * \brief get essential boundary conditions (only homogenous case is considered)
   * \ingroup maxwell_element
   * @return      error code
   */
  PetscErrorCode getEssentialBc() {
    MoABErrorCode rval;
    PetscFunctionBegin;
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,bit)) {
      if(bit->getName().compare(0,10,"ESSENTIALBC") == 0) {
        Range faces;
        rval = mField.get_moab().get_entities_by_type(bit->meshset,MBTRI,faces,true); CHKERRQ_MOAB(rval);
        rval = mField.get_moab().get_adjacencies(
          faces,1,true,blockData.essentialBc,moab::Interface::UNION
        ); CHKERRQ_MOAB(rval);
        blockData.essentialBc.merge(faces);
      }
    }
    if(blockData.essentialBc.empty()) {
      Range tets;
      rval = mField.get_moab().get_entities_by_type(0,MBTET,tets); CHKERRQ_MOAB(rval);
      Skinner skin(&mField.get_moab());
      Range skin_faces; // skin faces from 3d ents
      rval = skin.find_skin(0,tets,false,skin_faces); CHKERR_MOAB(rval);
      skin_faces = subtract(skin_faces,blockData.naturalBc);
      rval = mField.get_moab().get_adjacencies(
        skin_faces,1,true,blockData.essentialBc,moab::Interface::UNION
      ); CHKERRQ_MOAB(rval);
      blockData.essentialBc.merge(skin_faces);
    }
    PetscFunctionReturn(0);
  }

  /**
   * \brief build problem data structures
   * \ingroup maxwell_element
   * @return error code
   */
  PetscErrorCode createFields() {
    // MoABErrorCode rval;
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // Set entities bit level. each entity has bit level depending for example
    // on refinement level. In this case we do not refine mesh or not do
    // topological changes, simply set refinement level to zero on all entities.

    ierr = mField.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);

    // add fields
    ierr = mField.add_field(blockData.fieldName,HCURL,AINSWORTH_LOBBATO_BASE,1); CHKERRQ(ierr);
    ierr = mField.add_field("MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
    //meshset consisting all entities in mesh
    EntityHandle root_set = mField.get_moab().get_root_set();
    //add entities to field
    ierr = mField.add_ents_to_field_by_TETs(root_set,blockData.fieldName); CHKERRQ(ierr);

    // // The higher-order gradients can be gauged by locally skipping the
    // // corresponding degrees of freedom and basis functions in the higher-order
    // // edge-based, face-based and cell-based finite element subspaces.
    //
    // Range tris,edges;
    // rval = mField.get_moab().get_entities_by_type(root_set,MBTRI,tris,true); CHKERRQ_MOAB(rval);
    // rval = mField.get_moab().get_entities_by_type(root_set,MBEDGE,edges,true); CHKERRQ_MOAB(rval);
    //
    // // Set order in volume
    // Range bc_ents = unite(blockData.naturalBc,blockData.essentialBc);
    // Range vol_ents = subtract(unite(tris,edges),bc_ents);
    // ierr = mField.set_field_order(vol_ents,blockData.fieldName,blockData.oRder); CHKERRQ(ierr);
    // int gauged_order = 1;
    // ierr = mField.set_field_order(bc_ents,blockData.fieldName,gauged_order); CHKERRQ(ierr);

    // Set order on tets
    ierr = mField.set_field_order(root_set,MBTET,blockData.fieldName,blockData.oRder); CHKERRQ(ierr);
    ierr = mField.set_field_order(root_set,MBTRI,blockData.fieldName,blockData.oRder); CHKERRQ(ierr);
    ierr = mField.set_field_order(root_set,MBEDGE,blockData.fieldName,blockData.oRder); CHKERRQ(ierr);

    // Set geometry approximation ordered
    ierr = mField.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.set_field_order(root_set,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = mField.set_field_order(root_set,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = mField.set_field_order(root_set,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = mField.set_field_order(root_set,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

    //build field
    ierr = mField.build_fields(); CHKERRQ(ierr);

    // get HO geometry for 10 node tets
    // This method takes coordinates form edges mid nodes in 10 node tet and
    // project values on 2nd order hierarchical basis used to approx. geometry.
    Projection10NodeCoordsOnField ent_method_material(mField,"MESH_NODE_POSITIONS");
    ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  /**
   * \brief create finite elements
   * \ingroup maxwell_element
   *
   * Create volume and surface element. Surface element is used to integrate
   * natural boundary conditions.
   *
   * @return error code
   */
  PetscErrorCode createElements() {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    // //Elements
    ierr = mField.add_finite_element(blockData.feName); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(blockData.feName,blockData.fieldName); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(blockData.feName,blockData.fieldName); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(blockData.feName,blockData.fieldName); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(blockData.feName,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.add_finite_element(blockData.feNaturalBCName); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(blockData.feNaturalBCName,blockData.fieldName); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(blockData.feNaturalBCName,blockData.fieldName); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(blockData.feNaturalBCName,blockData.fieldName); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(blockData.feNaturalBCName,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TETs(0,blockData.feName); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(blockData.naturalBc,blockData.feNaturalBCName); CHKERRQ(ierr);
    //build finite elemnts
    ierr = mField.build_finite_elements(); CHKERRQ(ierr);
    //build adjacencies
    ierr = mField.build_adjacencies(BitRefLevel().set(0)); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
   * \brief create problem
   * \ingroup maxwell_element
   *
   * Problem is collection of finite elements. With the information on which
   * fields finite elements operates the matrix and  left and right hand side
   * vector could be created.
   *
   * Here we use Distributed mesh manager from PETSc as a inteface.
   *
   * @return error code
   */
  PetscErrorCode createProblem() {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    // set up DM
    DMType dm_name = "MAGNETIC_PROBLEM";
    ierr = DMRegister_MoFEM(dm_name); CHKERRQ(ierr);
    ierr = DMCreate(PETSC_COMM_WORLD,&blockData.dM);CHKERRQ(ierr);
    ierr = DMSetType(blockData.dM,dm_name);CHKERRQ(ierr);
    ierr = DMMoFEMCreateMoFEM(blockData.dM,&mField,dm_name,BitRefLevel().set(0)); CHKERRQ(ierr);
    ierr = DMSetFromOptions(blockData.dM); CHKERRQ(ierr);
    //add elements to blockData.dM
    ierr = DMMoFEMAddElement(blockData.dM,blockData.feName.c_str()); CHKERRQ(ierr);
    ierr = DMMoFEMAddElement(blockData.dM,blockData.feNaturalBCName.c_str()); CHKERRQ(ierr);
    ierr = DMSetUp(blockData.dM); CHKERRQ(ierr);
    // create matrices and vectors
    ierr = DMCreateGlobalVector(blockData.dM,&blockData.D); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
   * \brief destroy Distributed mesh manager
   * \ingroup maxwell_element
   * @return [description]
   */
  PetscErrorCode destroyProblem() {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = DMDestroy(&blockData.dM); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**  \brief solve problem
   * \ingroup maxwell_element
   *
   * Create matrices; integrate over elements; solve linear system of equations
   *
   */
  PetscErrorCode solveProblem() {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = DMCreateMatrix(blockData.dM,&blockData.A); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(blockData.dM,&blockData.F); CHKERRQ(ierr);
    ierr = VecDuplicate(blockData.F,&blockData.D); CHKERRQ(ierr);

    VolumeFE vol_fe(mField);
    vol_fe.getOpPtrVector().push_back(new OpCurlCurl(blockData));
    vol_fe.getOpPtrVector().push_back(new OpStab(blockData));
    TriFE tri_fe(mField);
    tri_fe.getOpPtrVector().push_back(new OpNaturalBC(blockData));

    ierr = MatZeroEntries(blockData.A); CHKERRQ(ierr);
    ierr = VecZeroEntries(blockData.F); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(blockData.F,ADD_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(blockData.F,ADD_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    ierr = DMoFEMLoopFiniteElements(blockData.dM,blockData.feName.c_str(),&vol_fe); CHKERRQ(ierr);
    ierr = DMoFEMLoopFiniteElements(blockData.dM,blockData.feNaturalBCName.c_str(),&tri_fe); CHKERRQ(ierr);

    ierr = MatAssemblyBegin(blockData.A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(blockData.A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(blockData.F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(blockData.F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(blockData.F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(blockData.F); CHKERRQ(ierr);

    // Boundary conditions
    std::vector<int> dofs_bc_indices;
    const MoFEM::MoFEMProblem *problem_ptr;
    ierr = DMMoFEMGetProblemPtr(blockData.dM,&problem_ptr); CHKERRQ(ierr);
    for(Range::iterator eit = blockData.essentialBc.begin();eit!=blockData.essentialBc.end();eit++) {
      for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(
        problem_ptr,blockData.fieldName,*eit,mField.getCommRank(),dof_ptr
      )) {
        std::bitset<8> pstatus(dof_ptr->get()->getPStatus());
        if(pstatus.test(0)) continue; //only local
        dofs_bc_indices.push_back(dof_ptr->get()->getPetscGlobalDofIdx());
      }
    }

    const double diag = 1;
    ierr = MatZeroRowsColumns(
      blockData.A,
      dofs_bc_indices.size(),
      dofs_bc_indices.empty()?PETSC_NULL:&*dofs_bc_indices.begin(),
      diag,
      PETSC_NULL,
      PETSC_NULL
    ); CHKERRQ(ierr);

    KSP solver;
    ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
    // ierr = KSPSetDM(solver,blockData.dM); CHKERRQ(ierr);
    ierr = KSPSetOperators(solver,blockData.A,blockData.A); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
    ierr = KSPSetUp(solver); CHKERRQ(ierr);
    ierr = KSPSolve(solver,blockData.F,blockData.D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(blockData.D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(blockData.D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = KSPDestroy(&solver); CHKERRQ(ierr);
    ierr = VecDestroy(&blockData.F); CHKERRQ(ierr);
    ierr = MatDestroy(&blockData.A); CHKERRQ(ierr);

    ierr = VecGhostUpdateBegin(blockData.D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(blockData.D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = DMoFEMMeshToLocalVector(
      blockData.dM,blockData.D,INSERT_VALUES,SCATTER_REVERSE
    ); CHKERRQ(ierr);
    ierr = VecDestroy(&blockData.D); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
   * \brief post-process results, i.e. save solution on the mesh
   * \ingroup maxwell_element
   * @return [description]
   */
  PetscErrorCode postProcessResults() {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    PostProcVolumeOnRefinedMesh post_proc(mField);
    ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc(blockData.fieldName); CHKERRQ(ierr);
    post_proc.getOpPtrVector().push_back(
      new OpPostProcessCurl(blockData,post_proc.postProcMesh,post_proc.mapGaussPts)
    );
    ierr = DMoFEMLoopFiniteElements(blockData.dM,blockData.feName.c_str(),&post_proc); CHKERRQ(ierr);
    ierr = post_proc.writeFile("out_values.h5m"); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


  /** \brief calculate and assemble CurlCurl matrix
   * \ingroup maxwell_element

  \f[
  \mathbf{A} = \int_\Omega \mu^{-1} \left( \nabla \times \mathbf{u}  \cdot \nabla \times \mathbf{v} \right) \textrm{d}\Omega
  \f]
  where
  \f[
  \mathbf{u} = \nabla \times \mathbf{B}
  \f]
  where \f$\mathbf{B}\f$ is magnetic flux and \f$\mu\f$ is magnetic permeability.

  For more details pleas look to \cite ivanyshyn2013computation

  */
  struct OpCurlCurl: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &blockData;
    OpCurlCurl(BlockData &data):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(data.fieldName,UserDataOperator::OPROWCOL),
    blockData(data) {
      sYmm = true;
    }

    MatrixDouble entityLocalMatrix;

    /**
     * \brief integrate matrix
     * @param  row_side local number of entity on element for row of the matrix
     * @param  col_side local number of entity on element for col of the matrix
     * @param  row_type type of row entity (EDGE/TRIANGLE/TETRAHEDRON)
     * @param  col_type type of col entity (EDGE/TRIANGLE/TETRAHEDRON)
     * @param  row_data structure of data, like base functions and associated methods to access those data on rows
     * @param  col_data structure of data, like base functions and associated methods to access those data on rows
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

      if(row_type==MBVERTEX) PetscFunctionReturn(0);
      if(col_type==MBVERTEX) PetscFunctionReturn(0);

      const int nb_row_dofs = row_data.getHcurlN().size2()/3;
      if(nb_row_dofs==0) PetscFunctionReturn(0);
      const int nb_col_dofs = col_data.getHcurlN().size2()/3;
      if(nb_col_dofs==0) PetscFunctionReturn(0);
      entityLocalMatrix.resize(nb_row_dofs,nb_col_dofs,false);
      entityLocalMatrix.clear();

      if(nb_row_dofs!=row_data.getFieldData().size()) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
          "Number of base functions and DOFs on entity is different on rows %d!=%d",
          nb_row_dofs,row_data.getFieldData().size()
        );
      }
      if(nb_col_dofs!=col_data.getFieldData().size()) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
          "Number of base functions and DOFs on entity is different on cols",
          nb_col_dofs,col_data.getFieldData().size()
        );
      }

      MatrixDouble row_curl_mat,col_curl_mat;
      FTensor::Index<'i',3> i;

      // cerr << row_data.getHcurlN() << endl;
      // cerr << row_data.getDiffHcurlN() << endl;

      const double c0 = 1./blockData.mU;
      const int nb_gauss_pts = row_data.getHcurlN().size1();

      for(int gg = 0;gg!=nb_gauss_pts;gg++) {

        // get integration weight scaled by volume
        double w = getGaussPts()(3,gg)*getVolume();
        // if ho geometry is given
        w *= getHoGaussPtsDetJac()(gg);
        ierr = getCurlOfHCurlBaseFunctions(
          row_side,row_type,row_data,gg,row_curl_mat
        ); CHKERRQ(ierr);
        ierr = getCurlOfHCurlBaseFunctions(
          col_side,col_type,col_data,gg,col_curl_mat
        ); CHKERRQ(ierr);

        // cerr << row_curl_mat << endl;
        // cerr << col_curl_mat << endl;

        FTensor::Tensor1<double*,3> t_row_curl(
          &row_curl_mat(0,HCURL0),&row_curl_mat(0,HCURL1),&row_curl_mat(0,HCURL2),3
        );
        for(int aa = 0;aa!=nb_row_dofs;aa++) {
          FTensor::Tensor0<double*> t_local_mat(&entityLocalMatrix(aa,0),1);
          FTensor::Tensor1<double*,3> t_col_curl(
            &col_curl_mat(0,HCURL0),&col_curl_mat(0,HCURL1),&col_curl_mat(0,HCURL2),3
          );
          for(int bb = 0;bb!=nb_col_dofs;bb++) {
            t_local_mat += c0*w*t_row_curl(i)*t_col_curl(i);
            ++t_col_curl;
            ++t_local_mat;
          }
          ++t_row_curl;
        }

      }

      // cerr << entityLocalMatrix << endl;
      // cerr << endl;

      ierr = MatSetValues(
        blockData.A,
        nb_row_dofs,&row_data.getIndices()[0],
        nb_col_dofs,&col_data.getIndices()[0],
        &entityLocalMatrix(0,0),ADD_VALUES
      ); CHKERRQ(ierr);

      if(row_side != col_side || row_type != col_type) {
        entityLocalMatrix = trans(entityLocalMatrix);
        ierr = MatSetValues(
          blockData.A,
          nb_col_dofs,&col_data.getIndices()[0],
          nb_row_dofs,&row_data.getIndices()[0],
          &entityLocalMatrix(0,0),ADD_VALUES
        ); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief calculate and assemble stabilization matrix
  * \ingroup maxwell_element

  \f[
  \mathbf{A} = \int_\Omega \epsilon \mathbf{u}  \cdot \mathbf{v} \textrm{d}\Omega
  \f]
  where \f$\epsilon\f$ is regularization parameter.

  */
  struct OpStab: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &blockData;
    OpStab(BlockData &data):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(data.fieldName,UserDataOperator::OPROWCOL),
    blockData(data) {
      sYmm = true;
    }

    MatrixDouble entityLocalMatrix;

    /**
     * \brief integrate matrix
     * @param  row_side local number of entity on element for row of the matrix
     * @param  col_side local number of entity on element for col of the matrix
     * @param  row_type type of row entity (EDGE/TRIANGLE/TETRAHEDRON)
     * @param  col_type type of col entity (EDGE/TRIANGLE/TETRAHEDRON)
     * @param  row_data structure of data, like base functions and associated methods to access those data on rows
     * @param  col_data structure of data, like base functions and associated methods to access those data on rows
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

      if(row_type==MBVERTEX) PetscFunctionReturn(0);
      if(col_type==MBVERTEX) PetscFunctionReturn(0);

      const int nb_row_dofs = row_data.getHcurlN().size2()/3;
      if(nb_row_dofs==0) PetscFunctionReturn(0);
      const int nb_col_dofs = col_data.getHcurlN().size2()/3;
      if(nb_col_dofs==0) PetscFunctionReturn(0);
      entityLocalMatrix.resize(nb_row_dofs,nb_col_dofs,false);
      entityLocalMatrix.clear();

      if(nb_row_dofs!=row_data.getFieldData().size()) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
          "Number of base functions and DOFs on entity is different on rows %d!=%d",
          nb_row_dofs,row_data.getFieldData().size()
        );
      }
      if(nb_col_dofs!=col_data.getFieldData().size()) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
          "Number of base functions and DOFs on entity is different on cols",
          nb_col_dofs,col_data.getFieldData().size()
        );
      }

      MatrixDouble row_curl_mat,col_curl_mat;
      FTensor::Index<'i',3> i;

      // cerr << row_data.getHcurlN() << endl;
      // cerr << row_data.getDiffHcurlN() << endl;

      const double c0 = 1./blockData.mU;
      const double c1 = blockData.ePsilon*c0;
      const int nb_gauss_pts = row_data.getHcurlN().size1();

      for(int gg = 0;gg!=nb_gauss_pts;gg++) {

        // get integration weight scaled by volume
        double w = getGaussPts()(3,gg)*getVolume();
        // if ho geometry is given
        w *= getHoGaussPtsDetJac()(gg);

        // cerr << row_curl_mat << endl;
        // cerr << col_curl_mat << endl;

        FTensor::Tensor1<const double*,3> t_row_base(
          &row_data.getHcurlN(gg)(0,HCURL0),
          &row_data.getHcurlN(gg)(0,HCURL1),
          &row_data.getHcurlN(gg)(0,HCURL2),3
        );

        for(int aa = 0;aa!=nb_row_dofs;aa++) {
          FTensor::Tensor0<double*> t_local_mat(&entityLocalMatrix(aa,0),1);
          FTensor::Tensor1<const double*,3> t_col_base(
            &col_data.getHcurlN(gg)(0,HCURL0),
            &col_data.getHcurlN(gg)(0,HCURL1),
            &col_data.getHcurlN(gg)(0,HCURL2),3
          );
          for(int bb = 0;bb!=nb_col_dofs;bb++) {
            t_local_mat += c1*w*t_row_base(i)*t_col_base(i);
            ++t_col_base;
            ++t_local_mat;
          }
          ++t_row_base;
        }

      }

      // cerr << entityLocalMatrix << endl;
      // cerr << endl;

      ierr = MatSetValues(
        blockData.A,
        nb_row_dofs,&row_data.getIndices()[0],
        nb_col_dofs,&col_data.getIndices()[0],
        &entityLocalMatrix(0,0),ADD_VALUES
      ); CHKERRQ(ierr);

      if(row_side != col_side || row_type != col_type) {
        entityLocalMatrix = trans(entityLocalMatrix);
        ierr = MatSetValues(
          blockData.A,
          nb_col_dofs,&col_data.getIndices()[0],
          nb_row_dofs,&row_data.getIndices()[0],
          &entityLocalMatrix(0,0),ADD_VALUES
        ); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

  };


  /** \brief calculate essential boundary conditions
    * \ingroup maxwell_element

    \f[
    \mathbf{A} = \int_{\partial\Omega} \ \mathbf{u}  \cdot \mathbf{j}_i \textrm{d}{\partial\Omega}
    \f]
    where \f$\mathbf{j}_i\f$ is current density function.

    Here simple current on coil is hard coded. In future more general implementation
    is needed.

    */
  struct OpNaturalBC: public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    BlockData &blockData;

    OpNaturalBC(BlockData &data):
    MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(data.fieldName,UserDataOperator::OPROW),
    blockData(data) {
    }

    VectorDouble naturalBC;

    /**
     * \brief integrate matrix
     * \ingroup maxwell_element
     * @param  row_side local number of entity on element for row of the matrix
     * @param  row_type type of row entity (EDGE/TRIANGLE/TETRAHEDRON)
     * @param  row_data structure of data, like base functions and associated methods to access those data on rows
     * @return          error code
     */
    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscErrorCode ierr;
      PetscFunctionBegin;

      if(row_type==MBVERTEX) PetscFunctionReturn(0);

      const int nb_row_dofs = row_data.getHcurlN().size2()/3;
      if(nb_row_dofs==0) PetscFunctionReturn(0);
      naturalBC.resize(nb_row_dofs,false);
      naturalBC.clear();

      FTensor::Index<'i',3> i;

      const int nb_gauss_pts = row_data.getHcurlN().size1();
      FTensor::Tensor1<double*,3> t_row_base = row_data.getFTensor1HcurlN<3>();

      FTensor::Tensor1<double*,3> t_tangent1 = getTensor1Tangent1AtGaussPt();
      FTensor::Tensor1<double*,3> t_tangent2 = getTensor1Tangent2AtGaussPt();

      for(int gg = 0;gg!=nb_gauss_pts;gg++) {

        // get integration weight scaled by volume
        double area;
        area = norm_2(getNormalsAtGaussPt(gg))*0.5;
        double w = getGaussPts()(2,gg)*area;

        // Current is on surface where natural bc are applied. It is set that
        // current is in XY plane, circular, around the coil.
        const double x = getHoCoordsAtGaussPts()(gg,0);
        const double y = getHoCoordsAtGaussPts()(gg,1);
        const double r = sqrt(x*x+y*y);
        FTensor::Tensor1<double,3> t_j;
        t_j(0) = -y/r;
        t_j(1) = +x/r;
        t_j(2) = 0;

        double a = t_j(i)*t_tangent1(i);
        double b = t_j(i)*t_tangent2(i);
        t_j(i) = a*t_tangent1(i)+b*t_tangent2(i);

        ++t_tangent1;
        ++t_tangent2;

        FTensor::Tensor0<double*> t_f(&naturalBC[0]);
        for(int aa = 0;aa!=nb_row_dofs;aa++) {
          t_f += w*t_row_base(i)*t_j(i);
          ++t_row_base;
          ++t_f;
        }

      }

      ierr = VecSetValues(
        blockData.F,
        row_data.getIndices().size(),
        &row_data.getIndices()[0],
        &naturalBC[0],
        ADD_VALUES
      ); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  /** \brief calculate and assemble CurlCurl matrix
    * \ingroup maxwell_element
    */
  struct OpPostProcessCurl: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &blockData;
    moab::Interface &postProcMesh;
    std::vector<EntityHandle> &mapGaussPts;

    OpPostProcessCurl(
      BlockData &data,
      moab::Interface &post_proc_mesh,
      std::vector<EntityHandle> &map_gauss_pts
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(
      data.fieldName,UserDataOperator::OPROW
    ),
    blockData(data),
    postProcMesh(post_proc_mesh),
    mapGaussPts(map_gauss_pts) {
    }
    virtual ~OpPostProcessCurl() {}

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscErrorCode ierr;
      MoABErrorCode rval;
      PetscFunctionBegin;

      if(row_type==MBVERTEX) PetscFunctionReturn(0);

      Tag th;
      double def_val[] = { 0,0,0 };
      rval = postProcMesh.tag_get_handle(
        "MAGNETIC_INDUCTION_FIELD",3,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
      ); CHKERRQ_MOAB(rval);
      const int nb_row_dofs = row_data.getHcurlN().size2()/3;
      if(nb_row_dofs==0) PetscFunctionReturn(0);
      const void* tags_ptr[mapGaussPts.size()];
      MatrixDouble row_curl_mat;
      FTensor::Index<'i',3> i;
      const int nb_gauss_pts = row_data.getHcurlN().size1();
      if(nb_gauss_pts!=mapGaussPts.size()) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
          "Inconsistency number of dofs %d!=%d",
          nb_gauss_pts,mapGaussPts.size()
        );
      }

      rval = postProcMesh.tag_get_by_ptr(th,&mapGaussPts[0],mapGaussPts.size(),tags_ptr); CHKERRQ_MOAB(rval);


      for(int gg = 0;gg!=nb_gauss_pts;gg++) {

        // get curl of base functions
        ierr = getCurlOfHCurlBaseFunctions(
          row_side,row_type,row_data,gg,row_curl_mat
        ); CHKERRQ(ierr);
        FTensor::Tensor1<double*,3> t_base_curl(
          &row_curl_mat(0,HCURL0),&row_curl_mat(0,HCURL1),&row_curl_mat(0,HCURL2),3
        );

        // get pointer to tag values on entity (i.e. vertex on refined post-processing mesh)
        double *ptr = &((double*)tags_ptr[gg])[0];
        FTensor::Tensor1<double*,3> t_curl(ptr,&ptr[1],&ptr[2]);

        // caclulate curl value
        for(int aa = 0;aa!=nb_row_dofs;aa++) {
          t_curl(i) += row_data.getFieldData()[aa]*t_base_curl(i);
          ++t_base_curl;
        }

      }
      PetscFunctionReturn(0);
    }

  };

};

#endif //__MAGNETICELEMENT_HPP__

/***************************************************************************//**
 * \defgroup maxwell_element Magnetic/Maxwell element
 * \ingroup user_modules
 ******************************************************************************/
