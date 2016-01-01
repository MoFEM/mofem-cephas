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

#include <MoFEM.hpp>
using namespace MoFEM;

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <adolc/adolc.h>
#include <NonLinearElasticElement.hpp>

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  Range CubitSIDESETs_meshsets;
  ierr = m_field.get_cubit_meshsets(SIDESET,CubitSIDESETs_meshsets); CHKERRQ(ierr);

  //ref meshset ref level 0
  ierr = m_field.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);
  //add entitities (by tets) to the field
  ierr = m_field.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);
  //set app. order
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }
  ierr = m_field.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);

  NonlinearElasticElement elastic(m_field,1);
  NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<double> double_kirchhoff_material;
  NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<adouble> adouble_kirchhoff_material;
  ierr = elastic.setBlocks(&double_kirchhoff_material,&adouble_kirchhoff_material); CHKERRQ(ierr);
  ierr = elastic.addElement("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = elastic.setOperators("SPATIAL_POSITION"); CHKERRQ(ierr);

  /*struct MyMat: public FunctionsToCalculatePiolaKirchhoffI {
    Interface& moAB;
    MyMat(Interface& moab): moAB(moab) {};
    PetscErrorCode calculateP_PiolaKirchhoffI(
     const BlockData block_data,
     const NumeredMoFEMFiniteElement *fe_ptr) {
     PetscFunctionBegin;
     //my stuff
     PetscFunctionReturn(0);
    }
  };
  MyMat mymat(moab);
  ierr = elastic.setOperators(mymat,"SPATIAL_POSITION"); CHKERRQ(ierr);*/

  //define problems
  ierr = m_field.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //set refinment level for problem
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
      if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
      EntityHandle ent = dof_ptr->get_ent();
      int dof_rank = dof_ptr->get_dof_coeff_idx();
      double &fval = dof_ptr->get_FieldData();
      if(node!=ent) {
        rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
        node = ent;
      }
      fval = scale_positions*coords[dof_rank];
    }
  }

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = m_field.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

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

  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"nonlinear_elastic.txt",&viewer); CHKERRQ(ierr);

  ierr = VecChop(F,1e-4); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecView(F,viewer); CHKERRQ(ierr);

  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);
  MatChop(Aij,1e-4);
  //MatView(Aij,PETSC_VIEWER_STDOUT_WORLD);
  MatView(Aij,viewer);
  //std::string wait;
  //std::cin >> wait;

  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = VecDestroy(&F); CHKERRQ(ierr);
  //ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);

  PetscFinalize();

  return 0;


}
