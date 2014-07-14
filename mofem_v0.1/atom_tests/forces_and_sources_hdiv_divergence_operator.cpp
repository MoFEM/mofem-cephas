#include "FEM.h"

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "ForcesAndSurcesCore.hpp"

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

#include<moab/Skinner.hpp>

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

  //create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;
  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //fields
  ierr = mField.add_field("HDIV",HDIV,1); CHKERRQ(ierr);
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"HDIV"); CHKERRQ(ierr);
  //set app. order
  int order = 1;
  ierr = mField.set_field_order(root_set,MBTET,"HDIV",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"HDIV",order); CHKERRQ(ierr);
  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //finite elements
  ierr = mField.add_finite_element("TET_FE"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("SKIN_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("TET_FE","HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("TET_FE","HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TET_FE","HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("SKIN_FE","HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("SKIN_FE","HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("SKIN_FE","HDIV"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = mField.add_ents_to_finite_element_by_TETs(root_set,"TET_FE"); CHKERRQ(ierr);
  Range tets;
  ierr = mField.get_entities_by_type_and_ref_level(BitRefLevel().set(0),BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
  Skinner skin(&moab);
  Range skin_faces; // skin faces from 3d ents
  rval = skin.find_skin(tets,false,skin_faces); CHKERR(rval);
  ierr = mField.add_ents_to_finite_element_by_TRIs(skin_faces,"SKIN_FE"); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finiteElementsPtr(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //problem
  ierr = mField.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("TEST_PROBLEM","TET_FE"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("TEST_PROBLEM","SKIN_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //mesh partitioning 
  //partition
  ierr = mField.simple_partition_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finiteElementsPtr("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  struct OpTetDivergence: public TetElementForcesAndSourcesCore::UserDataOperator {

    FieldData &dIv;
    OpTetDivergence(FieldData &div):
      TetElementForcesAndSourcesCore::UserDataOperator("HDIV"),dIv(div) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      if(data.getFieldData().size()==0) PetscFunctionReturn(0);

      cout << "type " << type << " side " << side << endl;

      int nb_gauss_pts = data.getDiffHdivN().size1();
      int nb_dofs = data.getFieldData().size();

      ublas::vector<FieldData> div_vec;
      div_vec.resize(nb_dofs,0);

      int gg = 0;
      for(;gg<nb_gauss_pts;gg++) {
	ierr = getDivergenceMatrixOperato_Hdiv(side,type,data,gg,div_vec); CHKERRQ(ierr);
	cout << div_vec << endl;
	cout << data.getDiffHdivN(gg) << endl;
      }

      cout << data.getDiffHdivN() << endl;

      cout << endl;

      PetscFunctionReturn(0);
    }	

  };

  struct MyFE: public TetElementForcesAndSourcesCore {

    MyFE(FieldInterface &m_field): TetElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return 0; }; //order/2; };

  };


  double divergence =0;

  MyFE tet_fe(mField);
  tet_fe.get_op_to_do_Rhs().push_back(new OpTetDivergence(divergence));

  ierr = mField.loop_finite_elements("TEST_PROBLEM","TET_FE",tet_fe);  CHKERRQ(ierr);
  
  cout << "divergence " << divergence << endl;


  ierr = PetscFinalize(); CHKERRQ(ierr);
}
