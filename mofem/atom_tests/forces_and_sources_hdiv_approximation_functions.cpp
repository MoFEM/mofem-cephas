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

//#include <DirichletBC.hpp>
//#include <PostProcOnRefMesh.hpp>

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //create one tet
  double tet_coords[] = {
    0,0,0,
    1,0,0,
    0,1,0,
    0,0,1
  };

  EntityHandle nodes[4];
  for(int nn = 0;nn<4;nn++) {
    rval = moab.create_vertex(&tet_coords[3*nn],nodes[nn]); CHKERRQ_MOAB(rval);
  }

  EntityHandle tet;
  rval = moab.create_element(MBTET,nodes,4,tet); CHKERRQ_MOAB(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("HDIV",HDIV,1); CHKERRQ(ierr);
  //ierr = m_field.add_field("L2",L2,1); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("TEST_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TEST_FE","HDIV"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TEST_FE","HDIV"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE","HDIV"); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);


  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"HDIV"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"TEST_FE"); CHKERRQ(ierr);

  //set app. order
  int order = 4;
  ierr = m_field.set_field_order(root_set,MBTET,"HDIV",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"HDIV",order); CHKERRQ(ierr);
  //ierr = m_field.set_field_order(root_set,MBTET,"L2",order); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning
  //partition
  ierr = m_field.partition_simple_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  typedef tee_device<std::ostream, std::ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;

  std::ofstream ofs("forces_and_sources_hdiv_approximation_functions.txt");
  TeeDevice my_tee(std::cout, ofs);
  TeeStream my_split(my_tee);


  struct OpPrintingHdivApproximationFunctions: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    TeeStream &mySplit;
    OpPrintingHdivApproximationFunctions(TeeStream &my_split):
    VolumeElementForcesAndSourcesCore::UserDataOperator("HDIV",UserDataOperator::OPROW),mySplit(my_split) {

    }

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data
    ) {
      PetscFunctionBegin;

      if(data.getFieldData().size()==0) PetscFunctionReturn(0);

      mySplit << std::endl << "type " << type << " side " << side << std::endl;
      mySplit.precision(5);

      const double eps = 1e-6;
      int dd = 0;
      int size = data.getHdivN().data().size();
      for(;dd<size;dd++) {
        if(fabs(data.getHdivN().data()[dd])<eps) data.getHdivN().data()[dd] = 0;
      }

      mySplit << std::fixed << data.getHdivN() << std::endl;


      PetscFunctionReturn(0);
    }

  };

  struct MyFE: public VolumeElementForcesAndSourcesCore {

    MyFE(MoFEM::Interface &m_field): VolumeElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return 1; };

  };

  MyFE tet_fe(m_field);
  tet_fe.getOpPtrVector().push_back(new OpPrintingHdivApproximationFunctions(my_split));

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE",tet_fe);  CHKERRQ(ierr);

  /*PostProcVolumeOnRefinedMesh post_proc(m_field);
  ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
  ierr = post_proc.addHdivFunctionsPostProc("HDIV");  CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE",post_proc);  CHKERRQ(ierr);
  rval = post_proc.postProcMesh.write_file("out.vtk","VTK",""); CHKERRQ_MOAB(rval);*/


  ierr = PetscFinalize(); CHKERRQ(ierr);
}
