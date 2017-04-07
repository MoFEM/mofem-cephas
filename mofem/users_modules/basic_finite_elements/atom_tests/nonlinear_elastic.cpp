/** \file nonlinear_elastic.cpp

 \brief Atom test for linear elastic dynamics.

 This is not exactly procedure for linear elastic dynamics, since jacobian is
 evaluated at every time step and snes procedure is involved. However it is
 implemented like that, to test methodology for general nonlinear problem.

*/

/* This file is part of MoFEM.
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

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

    moab::Core mb_instance;
    moab::Interface& moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
    ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
    }

    const char *option;
    option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    //ref meshset ref level 0
    ierr = m_field.seed_ref_level_3D(0,0); CHKERRQ(ierr);
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
    ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
    ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

    //Fields
    ierr = m_field.add_field("SPATIAL_POSITION",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
    //add entitities (by tets) to the field
    ierr = m_field.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);
    //set app. order
    PetscInt order;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      order = 1;
    }
    ierr = m_field.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);

    NonlinearElasticElement elastic(m_field,1);
    boost::shared_ptr<NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<double> >
    double_kirchhoff_material_ptr(new NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<double>());
    boost::shared_ptr<NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<adouble> >
    adouble_kirchhoff_material_ptr(new NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<adouble>());
    ierr = elastic.setBlocks(double_kirchhoff_material_ptr,adouble_kirchhoff_material_ptr); CHKERRQ(ierr);
    ierr = elastic.addElement("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = elastic.setOperators("SPATIAL_POSITION"); CHKERRQ(ierr);

    /*struct MyMat: public FunctionsToCalculatePiolaKirchhoffI {
    Interface& moAB;
    MyMat(Interface& moab): moAB(moab) {};
    PetscErrorCode calculateP_PiolaKirchhoffI(
    const BlockData block_data,
    const NumeredEntFiniteElement *fe_ptr) {
    PetscFunctionBegin;
    //my stuff
    PetscFunctionReturn(0);
  }
};
MyMat mymat(moab);
ierr = elastic.setOperators(mymat,"SPATIAL_POSITION"); CHKERRQ(ierr);*/

//define problems
ierr = m_field.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
//set refinement level for problem
ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);
//set finite elements for problems
ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);

//build field
ierr = m_field.build_fields(); CHKERRQ(ierr);

//use this to apply some strain field to the body (testing only)
double scale_positions = 2;
{
  EntityHandle node = 0;
  double coords[3];
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,"SPATIAL_POSITION",dof_ptr)) {
    if(dof_ptr->get()->getEntType()!=MBVERTEX) continue;
    EntityHandle ent = dof_ptr->get()->getEnt();
    int dof_rank = dof_ptr->get()->getDofCoeffIdx();
    double &fval = dof_ptr->get()->getFieldData();
    if(node!=ent) {
      rval = moab.get_coords(&ent,1,coords); CHKERRQ_MOAB(rval);
      node = ent;
    }
    fval = scale_positions*coords[dof_rank];
  }
}

//build finite elemnts
ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
//build adjacencies
ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);


ProblemsManager *prb_mng_ptr;
ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);
//build problem
ierr = prb_mng_ptr->buildProblem("ELASTIC_MECHANICS",true); CHKERRQ(ierr);
//partition
ierr = prb_mng_ptr->partitionProblem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
ierr = prb_mng_ptr->partitionFiniteElements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
ierr = prb_mng_ptr->partitionGhostDofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

//create matrices
Vec F;
ierr = m_field.VecCreateGhost("ELASTIC_MECHANICS",COL,&F); CHKERRQ(ierr);
//Vec D;
//ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
Mat Aij;
ierr = m_field.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

elastic.getLoopFeRhs().snes_f = F;
elastic.getLoopFeLhs().snes_B = Aij;

ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",elastic.getLoopFeRhs()); CHKERRQ(ierr);
ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",elastic.getLoopFeLhs()); CHKERRQ(ierr);
ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

// PetscViewer viewer;
// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"nonlinear_elastic.txt",&viewer); CHKERRQ(ierr);

// ierr = VecChop(F,1e-4); CHKERRQ(ierr);
// //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
// ierr = VecView(F,viewer); CHKERRQ(ierr);
//
// //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);
// MatChop(Aij,1e-4);
// //MatView(Aij,PETSC_VIEWER_STDOUT_WORLD);
// MatView(Aij,viewer);
// //std::string wait;
// //std::cin >> wait;
// ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

double sum = 0;
ierr = VecSum(F,&sum); CHKERRQ(ierr);
ierr = PetscPrintf(PETSC_COMM_WORLD,"sum  = %4.3e\n",sum); CHKERRQ(ierr);
double fnorm;
ierr = VecNorm(F,NORM_2,&fnorm); CHKERRQ(ierr);
ierr = PetscPrintf(PETSC_COMM_WORLD,"fnorm  = %9.8e\n",fnorm); CHKERRQ(ierr);

double mnorm;
ierr = MatNorm(Aij,NORM_1,&mnorm); CHKERRQ(ierr);
ierr = PetscPrintf(PETSC_COMM_WORLD,"mnorm  = %9.8e\n",mnorm); CHKERRQ(ierr);


if(fabs(sum)>1e-8) {
  SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Failed to pass test");
}
if(fabs(fnorm-5.12196914e+00)>1e-6) {
  SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Failed to pass test");
}
if(fabs(mnorm-5.48280139e+01)>1e-6) {
  SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Failed to pass test");
}


ierr = VecDestroy(&F); CHKERRQ(ierr);
//ierr = VecDestroy(&D); CHKERRQ(ierr);
ierr = MatDestroy(&Aij); CHKERRQ(ierr);


} catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFinalize();

  return 0;


}
