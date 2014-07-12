#include "FEM.h"

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "ForcesAndSurcesCore.hpp"
#include "PotsProcOnRefMesh.hpp"

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

  Core mb_instance;
  Interface& moab = mb_instance;
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
    rval = moab.create_vertex(&tet_coords[3*nn],nodes[nn]); CHKERR_PETSC(rval);
  }

  EntityHandle tet;
  rval = moab.create_element(MBTET,nodes,4,tet); CHKERR_PETSC(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("HDIV",HDIV,1); CHKERRQ(ierr);
  ierr = mField.add_field("L2",L2,1); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("TEST_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("TEST_FE","HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("TEST_FE","HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TEST_FE","HDIV"); CHKERRQ(ierr);

  //Problem
  ierr = mField.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);


  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"HDIV"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = mField.add_ents_to_finite_element_by_TETs(root_set,"TEST_FE"); CHKERRQ(ierr);

  //set app. order
  int order = 4;
  ierr = mField.set_field_order(root_set,MBTET,"HDIV",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"HDIV",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTET,"L2",order); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = mField.build_finiteElementsPtr(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning 
  //partition
  ierr = mField.simple_partition_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finiteElementsPtr("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  typedef tee_device<ostream, ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;

  ofstream ofs("forces_and_sources_hdiv_approximation_functions.txt");
  TeeDevice my_tee(cout, ofs); 
  TeeStream my_split(my_tee);


  struct OpPrintingHdivApproximationFunctions: public TetElementForcesAndSourcesCore::UserDataOperator {

    TeeStream &mySplit;
    OpPrintingHdivApproximationFunctions(TeeStream &my_split):
      TetElementForcesAndSourcesCore::UserDataOperator("HDIV"),mySplit(my_split) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getFieldData().size()==0) PetscFunctionReturn(0);

      mySplit << endl << "type " << type << " side " << side << endl;
      mySplit.precision(5);
      mySplit << std::fixed << data.getDiffN() << endl;
      

      PetscFunctionReturn(0);
    }	

  };

  struct MyFE: public TetElementForcesAndSourcesCore {

    MyFE(FieldInterface &m_field): TetElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return 1; };

  };

  MyFE tet_fe(mField);
  tet_fe.get_op_to_do_Rhs().push_back(new OpPrintingHdivApproximationFunctions(my_split));

  ierr = mField.loop_finite_elements("TEST_PROBLEM","TEST_FE",tet_fe);  CHKERRQ(ierr);

  PostPocOnRefinedMesh post_proc(mField);
  ierr = post_proc.generateRefereneElemenMesh(); CHKERRQ(ierr);
  ierr = post_proc.addHdivFunctionsPostProc("HDIV");  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("TEST_PROBLEM","TEST_FE",post_proc);  CHKERRQ(ierr);

  rval = post_proc.postProcMesh.write_file("out.vtk","VTK",""); CHKERR_PETSC(rval);


  ierr = PetscFinalize(); CHKERRQ(ierr);
}
