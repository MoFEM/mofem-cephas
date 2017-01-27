/** \file stability.cpp
 * \ingroup nonlinear_elastic_elem
 *
 * Solves stability problem. Currently uses 3d tetrahedral elements.
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

#include <BasicFiniteElements.hpp>
using namespace MoFEM;
#include <Hooke.hpp>

#undef EPS
#include <slepceps.h>

#include <SurfacePressureComplexForLazy.hpp>

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

template<typename TYPE>
struct MyMat_double: public NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<TYPE> {

  bool doAotherwiseB;
  MyMat_double(): doAotherwiseB(true) {};

  ublas::matrix<double> D_lambda,D_mu,D;
  ublas::vector<TYPE> sTrain,sTrain0,sTress;
  ublas::matrix<adouble> invF,CauchyStress;

  virtual PetscErrorCode calculateP_PiolaKirchhoffI(
    const NonlinearElasticElement::BlockData block_data,
    const NumeredEntFiniteElement *fe_ptr) {
    PetscFunctionBegin;

    try {

      double lambda = LAMBDA(block_data.E,block_data.PoissonRatio);
      double mu = MU(block_data.E,block_data.PoissonRatio);
      if(D_lambda.size1()==0) {
        D_lambda.resize(6,6);
        D_lambda.clear();
        for(int rr = 0;rr<3;rr++) {
          for(int cc = 0;cc<3;cc++) {
            D_lambda(rr,cc) = 1;
          }
        }
      }
      if(D_mu.size1()==0) {
        D_mu.resize(6,6);
        D_mu.clear();
        for(int rr = 0;rr<6;rr++) {
          D_mu(rr,rr) = rr<3 ? 2 : 1;
        }
      }
      D.resize(6,6);
      noalias(D) = lambda*D_lambda + mu*D_mu;

      if(doAotherwiseB) {
        sTrain.resize(6);
        sTrain[0] = this->F(0,0)-1;
        sTrain[1] = this->F(1,1)-1;
        sTrain[2] = this->F(2,2)-1;
        sTrain[3] = this->F(0,1)+this->F(1,0);
        sTrain[4] = this->F(1,2)+this->F(2,1);
        sTrain[5] = this->F(0,2)+this->F(2,0);
        sTress.resize(6);
        noalias(sTress) = prod(D,sTrain);
        this->P.resize(3,3);
        this->P(0,0) = sTress[0];
        this->P(1,1) = sTress[1];
        this->P(2,2) = sTress[2];
        this->P(0,1) = this->P(1,0) = sTress[3];
        this->P(1,2) = this->P(2,1) = sTress[4];
        this->P(0,2) = this->P(2,0) = sTress[5];
        //std::cerr << this->P << std::endl;
      } else {
        adouble J;
        ierr = this->dEterminatnt(this->F,J); CHKERRQ(ierr);
        invF.resize(3,3);
        ierr = this->iNvert(J,this->F,invF); CHKERRQ(ierr);
        sTrain0.resize(6,0);
        noalias(sTress) = prod(D,sTrain0);
        CauchyStress.resize(3,3);
        CauchyStress(0,0) = sTress[0];
        CauchyStress(1,1) = sTress[1];
        CauchyStress(2,2) = sTress[2];
        CauchyStress(0,1) = CauchyStress(1,0) = sTress[3];
        CauchyStress(1,2) = CauchyStress(2,1) = sTress[4];
        CauchyStress(0,2) = CauchyStress(2,0) = sTress[5];
        //std::cerr << D << std::endl;
        //std::cerr << CauchyStress << std::endl;
        noalias(this->P) = J*prod(CauchyStress,trans(invF));
      }

    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "throw in method: " << ex.what() << std::endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }

};

template<typename TYPE>
struct MyMat: public MyMat_double<TYPE> {

  int nbActiveVariables0;

  virtual PetscErrorCode setUserActiveVariables(
    int &nb_active_variables) {
    PetscFunctionBegin;

    try {

      this->sTrain0.resize(6);
      ublas::matrix<double> &G0 = (this->commonDataPtr->gradAtGaussPts["D0"][this->gG]);
      this->sTrain0[0] <<= G0(0,0);
      this->sTrain0[1] <<= G0(1,1);
      this->sTrain0[2] <<= G0(2,2);
      this->sTrain0[3] <<= (G0(1,0) + G0(0,1));
      this->sTrain0[4] <<= (G0(2,1) + G0(1,2));
      this->sTrain0[5] <<= (G0(2,0) + G0(0,2));
      nbActiveVariables0 = nb_active_variables;
      nb_active_variables += 6;

    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "throw in method: " << ex.what() << std::endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode setUserActiveVariables(
    ublas::vector<double> &active_varibles) {
    PetscFunctionBegin;

    try {

      int shift = nbActiveVariables0; // is a number of elements in F
      ublas::matrix<double> &G0 = (this->commonDataPtr->gradAtGaussPts["D0"][this->gG]);
      active_varibles[shift+0] = G0(0,0);
      active_varibles[shift+1] = G0(1,1);
      active_varibles[shift+2] = G0(2,2);
      active_varibles[shift+3] = G0(0,1)+G0(1,0);
      active_varibles[shift+4] = G0(1,2)+G0(2,1);
      active_varibles[shift+5] = G0(0,2)+G0(2,0);

    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "throw in method: " << ex.what() << std::endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }


};

int main(int argc, char *argv[]) {

  //PetscInitialize(&argc,&argv,(char *)0,help);
  SlepcInitialize(&argc,&argv,(char*)0,help);

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);


  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  // use this if your mesh is partotioned and you run code on parts,
  // you can solve very big problems
  PetscBool is_partitioned = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-my_is_partitioned",&is_partitioned,&flg); CHKERRQ(ierr);

  if(is_partitioned == PETSC_TRUE) {
    //Read mesh to MOAB
    const char *option;
    option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
    rval = pcomm->resolve_shared_ents(0,3,0); CHKERRQ_MOAB(rval);
    rval = pcomm->resolve_shared_ents(0,3,1); CHKERRQ_MOAB(rval);
    rval = pcomm->resolve_shared_ents(0,3,2); CHKERRQ_MOAB(rval);
  } else {
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  }

  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  Range CubitSIDESETs_meshsets;
  ierr = m_field.get_cubit_meshsets(SIDESET,CubitSIDESETs_meshsets); CHKERRQ(ierr);

  //ref meshset ref level 0
  ierr = m_field.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  bool check_if_spatial_field_exist = m_field.check_field("SPATIAL_POSITION");
  ierr = m_field.add_field("SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("EIGEN_VECTOR",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("D0",H1,3,MF_ZERO); CHKERRQ(ierr);

  //add entitities (by tets) to the field
  ierr = m_field.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"EIGEN_VECTOR"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"D0"); CHKERRQ(ierr);

  Hooke<double> mat_double;
  MyMat<adouble> mat_adouble;

  NonlinearElasticElement elastic(m_field,2);
  ierr = elastic.setBlocks(&mat_double,&mat_adouble); CHKERRQ(ierr);
  ierr = elastic.addElement("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC","EIGEN_VECTOR"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC","D0"); CHKERRQ(ierr);

  elastic.feRhs.getOpPtrVector().push_back(
    new NonlinearElasticElement::OpGetCommonDataAtGaussPts("D0",elastic.commonData));
  elastic.feLhs.getOpPtrVector().push_back(
    new NonlinearElasticElement::OpGetCommonDataAtGaussPts("D0",elastic.commonData));
  ierr = elastic.setOperators("SPATIAL_POSITION"); CHKERRQ(ierr);

  //define problems
  ierr = m_field.add_problem("ELASTIC_MECHANICS",MF_ZERO); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  //set app. order

  PetscInt disp_order;
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&disp_order,&flg); CHKERRQ(ierr);
  if(flg!=PETSC_TRUE) {
    disp_order = 1;
  }

  ierr = m_field.set_field_order(0,MBTET,"SPATIAL_POSITION",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_POSITION",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_POSITION",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"EIGEN_VECTOR",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"EIGEN_VECTOR",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"EIGEN_VECTOR",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"EIGEN_VECTOR",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"D0",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"D0",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"D0",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"D0",1); CHKERRQ(ierr);

  ierr = m_field.add_finite_element("NEUAMNN_FE",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("NEUAMNN_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","NEUAMNN_FE"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    Range tris;
    rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERRQ_MOAB(rval);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(tris,"NEUAMNN_FE"); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
    Range tris;
    rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERRQ_MOAB(rval);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(tris,"NEUAMNN_FE"); CHKERRQ(ierr);
  }
  //add nodal force element
  ierr = MetaNodalForces::addElement(m_field,"SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","FORCE_FE"); CHKERRQ(ierr);

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //10 node tets
  if(!check_if_spatial_field_exist) {
    Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
    ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
    Projection10NodeCoordsOnField ent_method_spatial(m_field,"SPATIAL_POSITION");
    ierr = m_field.loop_dofs("SPATIAL_POSITION",ent_method_spatial); CHKERRQ(ierr);
    //ierr = m_field.set_field(0,MBTRI,"SPATIAL_POSITION"); CHKERRQ(ierr);
    //ierr = m_field.set_field(0,MBTET,"SPATIAL_POSITION"); CHKERRQ(ierr);
    //ierr = m_field.field_axpy(1,"SPATIAL_POSITION","D0",true); CHKERRQ(ierr);
  }

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build database
  ProblemsManager *prb_mng_ptr;
  ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);
  if(is_partitioned) {
    ierr = prb_mng_ptr->buildProblemOnDistributedMesh(true); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements("ELASTIC_MECHANICS",true,0,pcomm->size(),1); CHKERRQ(ierr);
  } else {
    ierr = prb_mng_ptr->buildProblem(); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionProblem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  }
  ierr = prb_mng_ptr->partitionGhostDofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //create matrices
  Vec F;
  ierr = m_field.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F); CHKERRQ(ierr);
  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  Mat Aij;
  ierr = m_field.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  //surface forces
  NeummanForcesSurfaceComplexForLazy neumann_forces(m_field,Aij,F);
  NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &neumann = neumann_forces.getLoopSpatialFe();
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = neumann.addForce(it->getMeshsetId()); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
    ierr = neumann.addPreassure(it->getMeshsetId()); CHKERRQ(ierr);
  }
  SpatialPositionsBCFEMethodPreAndPostProc my_Dirichlet_bc(m_field,"SPATIAL_POSITION",Aij,D,F);

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);

  ierr = m_field.set_local_ghost_vector("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //F Vector
  //preproc
  my_Dirichlet_bc.snes_ctx = SnesMethod::CTX_SNESSETFUNCTION;
  my_Dirichlet_bc.snes_x = D;
  my_Dirichlet_bc.snes_f = F;
  ierr = m_field.problem_basic_method_preProcess("ELASTIC_MECHANICS",my_Dirichlet_bc); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = m_field.set_local_ghost_vector("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //elem loops
  //noadl forces
  boost::ptr_map<std::string,NodalForce> nodal_forces;
  string fe_name_str ="FORCE_FE";
  nodal_forces.insert(fe_name_str,new NodalForce(m_field));
  ierr = MetaNodalForces::setOperators(m_field,nodal_forces,F,"SPATIAL_POSITION"); CHKERRQ(ierr);
  boost::ptr_map<std::string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS",fit->first,fit->second->getLoopFe()); CHKERRQ(ierr);
  }
  //surface forces
  neumann.snes_ctx = SnesMethod::CTX_SNESSETFUNCTION;
  neumann.snes_x = D;
  neumann.snes_f = F;
  m_field.loop_finite_elements("ELASTIC_MECHANICS","NEUAMNN_FE",neumann);
  //stiffnes
  elastic.getLoopFeRhs().snes_ctx = SnesMethod::CTX_SNESSETFUNCTION;
  elastic.getLoopFeRhs().snes_x = D;
  elastic.getLoopFeRhs().snes_f = F;
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",elastic.getLoopFeRhs()); CHKERRQ(ierr);
  //postproc
  ierr = m_field.problem_basic_method_postProcess("ELASTIC_MECHANICS",my_Dirichlet_bc); CHKERRQ(ierr);

  //Aij Matrix
  //preproc
  my_Dirichlet_bc.snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
  my_Dirichlet_bc.snes_B = Aij;
  ierr = m_field.problem_basic_method_preProcess("ELASTIC_MECHANICS",my_Dirichlet_bc); CHKERRQ(ierr);
  //surface forces
  //neumann.snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
  //neumann.snes_B = Aij;
  //ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","NEUAMNN_FE",neumann); CHKERRQ(ierr);
  //stiffnes
  elastic.getLoopFeLhs().snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
  elastic.getLoopFeLhs().snes_B = Aij;
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",elastic.getLoopFeLhs()); CHKERRQ(ierr);
  //postproc
  ierr = m_field.problem_basic_method_postProcess("ELASTIC_MECHANICS",my_Dirichlet_bc); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  //Matrix View
  //MatView(Aij,PETSC_VIEWER_STDOUT_WORLD);
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);

  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = VecZeroEntries(D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = m_field.set_other_global_ghost_vector(
    "ELASTIC_MECHANICS","SPATIAL_POSITION","D0",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  Mat Bij;
  ierr = MatDuplicate(Aij,MAT_SHARE_NONZERO_PATTERN,&Bij); CHKERRQ(ierr);
  //ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  ierr = MatZeroEntries(Bij); CHKERRQ(ierr);

  /*//Aij Matrix
  //preproc
  my_Dirichlet_bc.snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
  my_Dirichlet_bc.snes_B = Aij;
  ierr = m_field.problem_basic_method_preProcess("ELASTIC_MECHANICS",my_Dirichlet_bc); CHKERRQ(ierr);
  //stiffnes
  elastic.getLoopFeLhs().snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
  elastic.getLoopFeLhs().snes_B = Aij;
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",elastic.getLoopFeLhs()); CHKERRQ(ierr);
  //postproc
  ierr = m_field.problem_basic_method_postProcess("ELASTIC_MECHANICS",my_Dirichlet_bc); CHKERRQ(ierr);*/

  //Bij Matrix
  mat_adouble.doAotherwiseB = false;
  //preproc
  my_Dirichlet_bc.snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
  my_Dirichlet_bc.snes_B = Bij;
  ierr = m_field.problem_basic_method_preProcess("ELASTIC_MECHANICS",my_Dirichlet_bc); CHKERRQ(ierr);
  //surface forces
  neumann.snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
  neumann.snes_B = Bij;
  PetscBool is_conservative = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-my_is_conservative",&is_conservative,&flg); CHKERRQ(ierr);
  if(is_conservative) {
    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","NEUAMNN_FE",neumann); CHKERRQ(ierr);
  }
  //stiffnes
  elastic.getLoopFeLhs().snes_ctx = SnesMethod::CTX_SNESSETJACOBIAN;
  elastic.getLoopFeLhs().snes_B = Bij;
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",elastic.getLoopFeLhs()); CHKERRQ(ierr);
  //postproc
  ierr = m_field.problem_basic_method_postProcess("ELASTIC_MECHANICS",my_Dirichlet_bc); CHKERRQ(ierr);

  ierr = MatSetOption(Bij,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Bij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Bij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  //Matrix View
  //MatView(Bij,PETSC_VIEWER_STDOUT_WORLD);
  //MatView(Bij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;

  EPS eps;
  ST st;
  EPSType type;
  PetscReal tol;
  PetscInt nev,maxit,its,lits;

  /*
    Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps); CHKERRQ(ierr);
  /*
    Set operators. In this case, it is a generalized eigenvalue problem
  */
  ierr = EPSSetOperators(eps,Bij,Aij); CHKERRQ(ierr);
  /*
    Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);
  /*
    Optional: Get some information from the solver and display it
  */
  ierr = EPSSolve(eps); CHKERRQ(ierr);

  /*
    Optional: Get some information from the solver and display it
  */
  ierr = EPSGetIterationNumber(eps,&its); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);
  ierr = EPSGetST(eps,&st); CHKERRQ(ierr);
  ierr = STGetOperationCounters(st,NULL,&lits); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %D\n",lits);
  ierr = EPSGetType(eps,&type); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n",type);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);
  ierr = EPSGetTolerances(eps,&tol,&maxit); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);

  //get solutions
  PostProcVolumeOnRefinedMesh post_proc(m_field);
  ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("EIGEN_VECTOR"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("EIGEN_VECTOR"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("D0"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("D0"); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",post_proc); CHKERRQ(ierr);
  ierr = post_proc.writeFile("out.h5m"); CHKERRQ(ierr);

  PetscScalar eigr,eigi,nrm2r;
  for(int nn = 0;nn<nev;nn++) {
    ierr = EPSGetEigenpair(eps,nn,&eigr,&eigi,D,PETSC_NULL); CHKERRQ(ierr);
    ierr = VecNorm(D,NORM_2,&nrm2r); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD," ncov = %D eigr = %.4g eigi = %.4g (inv eigr = %.4g) nrm2r = %.4g\n",nn,eigr,eigi,1./eigr,nrm2r);
    std::ostringstream o1;
    o1 << "eig_" << nn << ".h5m";
    ierr = m_field.set_other_global_ghost_vector(
      "ELASTIC_MECHANICS","SPATIAL_POSITION","EIGEN_VECTOR",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",post_proc); CHKERRQ(ierr);
    ierr = post_proc.writeFile(o1.str().c_str()); CHKERRQ(ierr);
  }

  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = MatDestroy(&Bij); CHKERRQ(ierr);
  ierr = EPSDestroy(&eps); CHKERRQ(ierr);

  SlepcFinalize();
  //PetscFinalize();

  return 0;
}
