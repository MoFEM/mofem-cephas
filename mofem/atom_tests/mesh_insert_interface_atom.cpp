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

static char help[] = "teting interface inserting algorithm\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc,&argv,(char *)0,help);

  try {

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRG(ierr);
  #else
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRG(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  /*char mesh_out_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_out_file",mesh_out_file_name,255,&flg); CHKERRG(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_out_file (MESH FILE NEEDED)");
  }**/

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  const char *option;
  option = "";//"PARALLEL=BCAST";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  PrismInterface *interface;
  ierr = m_field.getInterface(interface); CHKERRG(ierr);

  ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(0,3,BitRefLevel().set(0)); CHKERRG(ierr);
  std::vector<BitRefLevel> bit_levels;
  bit_levels.push_back(BitRefLevel().set(0));

  int ll = 1;
  //for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|INTERFACESET,cit)) {
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,cit)) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Insert Interface %d\n",cit->getMeshsetId()); CHKERRG(ierr);
    EntityHandle cubit_meshset = cit->getMeshset();
    {
      //get tet enties form back bit_level
      EntityHandle ref_level_meshset = 0;
      rval = moab.create_meshset(MESHSET_SET,ref_level_meshset); CHKERRG(rval);
      ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(bit_levels.back(),BitRefLevel().set(),MBTET,ref_level_meshset); CHKERRG(ierr);
      ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(bit_levels.back(),BitRefLevel().set(),MBPRISM,ref_level_meshset); CHKERRG(ierr);
      Range ref_level_tets;
      rval = moab.get_entities_by_handle(ref_level_meshset,ref_level_tets,true); CHKERRG(rval);
      //get faces and test to split
      ierr = interface->getSides(cubit_meshset,bit_levels.back(),true,0); CHKERRG(ierr);
      //set new bit level
      bit_levels.push_back(BitRefLevel().set(ll++));
      //split faces and
      ierr = interface->splitSides(ref_level_meshset,bit_levels.back(),cubit_meshset,true,true,0); CHKERRG(ierr);
      //clean meshsets
      rval = moab.delete_entities(&ref_level_meshset,1); CHKERRG(rval);
    }
    //update cubit meshsets
    for(_IT_CUBITMESHSETS_FOR_LOOP_(m_field,ciit)) {
      EntityHandle cubit_meshset = ciit->meshset;
      ierr = m_field.getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(cubit_meshset,bit_levels.back(),cubit_meshset,MBVERTEX,true); CHKERRG(ierr);
      ierr = m_field.getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(cubit_meshset,bit_levels.back(),cubit_meshset,MBEDGE,true); CHKERRG(ierr);
      ierr = m_field.getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(cubit_meshset,bit_levels.back(),cubit_meshset,MBTRI,true); CHKERRG(ierr);
      ierr = m_field.getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(cubit_meshset,bit_levels.back(),cubit_meshset,MBTET,true); CHKERRG(ierr);
    }
  }

  //add filds
  ierr = m_field.add_field("H1FIELD_SCALAR",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRG(ierr);

  //add finite elements
  ierr = m_field.add_finite_element("ELEM_SCALAR"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_row("ELEM_SCALAR","H1FIELD_SCALAR"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_col("ELEM_SCALAR","H1FIELD_SCALAR"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELEM_SCALAR","H1FIELD_SCALAR"); CHKERRG(ierr);
  //FE Interface
  ierr = m_field.add_finite_element("INTERFACE"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_row("INTERFACE","H1FIELD_SCALAR"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_col("INTERFACE","H1FIELD_SCALAR"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("INTERFACE","H1FIELD_SCALAR"); CHKERRG(ierr);

  //add ents to field and set app. order
  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"H1FIELD_SCALAR"); CHKERRG(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"H1FIELD_SCALAR",1); CHKERRG(ierr);

  //add finite elements entities
  //all TETS and PRIMS are added to finite elements, for testin pruposes.
  //in some practiacl applications to save memory, you would like to add elements
  //from particular refinement level (see: m_field.add_ents_to_finite_element_by_bit_ref(...)
  ierr = m_field.add_ents_to_finite_element_by_type(0,MBTET,"ELEM_SCALAR",true); CHKERRG(ierr);
  ierr = m_field.add_ents_to_finite_element_by_type(0,MBPRISM,"INTERFACE",true); CHKERRG(ierr);

  //add problems
  //set problem for all last two levels, only for testing pruposes
  for(int lll = ll-2;lll<ll;lll++) {
    std::stringstream problem_name;
    problem_name << "PROBLEM_SCALAR_" << lll;
    ierr = m_field.add_problem(problem_name.str()); CHKERRG(ierr);
    //define problems and finite elements
    ierr = m_field.modify_problem_add_finite_element(problem_name.str(),"ELEM_SCALAR"); CHKERRG(ierr);
    ierr = m_field.modify_problem_add_finite_element(problem_name.str(),"INTERFACE"); CHKERRG(ierr);
  }

  //set problem level
  for(int lll = ll-2;lll<ll;lll++) {
    std::stringstream problem_name;
    problem_name << "PROBLEM_SCALAR_" << lll;
    std::stringstream message;
    message << "set problem problem < " << problem_name.str() << " > bit level " << bit_levels[lll] << std::endl;
    ierr = PetscPrintf(PETSC_COMM_WORLD,message.str().c_str()); CHKERRG(ierr);
    ierr = m_field.modify_problem_ref_level_add_bit(problem_name.str(),bit_levels[lll]); CHKERRG(ierr);
  }

  //build fields
  ierr = m_field.build_fields(); CHKERRG(ierr);
  //build finite elements
  ierr = m_field.build_finite_elements(); CHKERRG(ierr);
  //build adjacencies
  //Its build adjacencies for all ements in databse,
  //for pratical applications consider to build adjacencies
  //only for refinemnt levels which you use for calulations
  ierr = m_field.build_adjacencies(BitRefLevel().set()); CHKERRG(ierr);

  Range tets_back_bit_level;
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(bit_levels.back(),BitRefLevel().set(),tets_back_bit_level); CHKERRG(ierr);

  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,cit)) {

    EntityHandle cubit_meshset = cit->getMeshset();

    BlockSetAttributes mydata;
    ierr = cit->getAttributeDataStructure(mydata); CHKERRG(ierr);
    std::cout << mydata << std::endl;

    Range tets;
    rval = moab.get_entities_by_type(cubit_meshset,MBTET,tets,true); CHKERRG(rval);
    tets = intersect(tets_back_bit_level,tets);
    Range nodes;
    rval = moab.get_connectivity(tets,nodes,true); CHKERRG(rval);

    for(Range::iterator nit = nodes.begin(); nit!=nodes.end(); nit++) {
      double coords[3];
      rval = moab.get_coords(&*nit,1,coords); CHKERRG(rval);
      coords[0] += mydata.data.User1;
      coords[1] += mydata.data.User2;
      coords[2] += mydata.data.User3;
      rval = moab.set_coords(&*nit,1,coords); CHKERRG(rval);
    }

  }


  ProblemsManager *prb_mng_ptr;
  ierr = m_field.getInterface(prb_mng_ptr); CHKERRG(ierr);

  //partition
  for(int lll = ll-2;lll<ll;lll++) {
    //build problem
    std::stringstream problem_name;
    problem_name << "PROBLEM_SCALAR_" << lll;
    ierr = prb_mng_ptr->buildProblem(problem_name.str(),true); CHKERRG(ierr);
    ierr = prb_mng_ptr->partitionProblem(problem_name.str()); CHKERRG(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements(problem_name.str()); CHKERRG(ierr);
    ierr = prb_mng_ptr->partitionGhostDofs(problem_name.str()); CHKERRG(ierr);
  }

  std::ofstream myfile;
  myfile.open("mesh_insert_interface.txt");

  EntityHandle out_meshset_tet;
  rval = moab.create_meshset(MESHSET_SET,out_meshset_tet); CHKERRG(rval);


  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(bit_levels.back(),BitRefLevel().set(),MBTET,out_meshset_tet); CHKERRG(ierr);
  Range tets;
  rval = moab.get_entities_by_handle(out_meshset_tet,tets,true); CHKERRG(rval);
  for(Range::iterator tit = tets.begin();tit!=tets.end();tit++) {
    int num_nodes;
    const EntityHandle* conn;
    rval = moab.get_connectivity(*tit,conn,num_nodes,true); CHKERRG(rval);

    for(int nn = 0;nn<num_nodes;nn++) {
      std::cout << conn[nn] << " ";
      myfile << conn[nn] << " ";
    }
    std::cout << std::endl;
    myfile << std::endl;

  }
  EntityHandle out_meshset_prism;
  rval = moab.create_meshset(MESHSET_SET,out_meshset_prism); CHKERRG(rval);
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(bit_levels.back(),BitRefLevel().set(),MBPRISM,out_meshset_prism); CHKERRG(ierr);
  Range prisms;
  rval = moab.get_entities_by_handle(out_meshset_prism,prisms); CHKERRG(rval);
  for(Range::iterator pit = prisms.begin();pit!=prisms.end();pit++) {
    int num_nodes;
    const EntityHandle* conn;
    rval = moab.get_connectivity(*pit,conn,num_nodes,true); CHKERRG(rval);

    for(int nn = 0;nn<num_nodes;nn++) {
      std::cout << conn[nn] << " ";
      myfile << conn[nn] << " ";
    }
    std::cout << std::endl;
    myfile << std::endl;

  }
  myfile.close();

  rval = moab.write_file("out_tet.vtk","VTK","",&out_meshset_tet,1); CHKERRG(rval);
  rval = moab.write_file("out_prism.vtk","VTK","",&out_meshset_prism,1); CHKERRG(rval);

  EntityHandle out_meshset_tets_and_prism;
  rval = moab.create_meshset(MESHSET_SET,out_meshset_tets_and_prism); CHKERRG(rval);
  rval = moab.add_entities(out_meshset_tets_and_prism,tets); CHKERRG(rval);
  rval = moab.add_entities(out_meshset_tets_and_prism,prisms); CHKERRG(rval);
  rval = moab.write_file("out_tets_and_prisms.vtk","VTK","",&out_meshset_tets_and_prism,1); CHKERRG(rval);

  EntityHandle out_meshset_tris;
  rval = moab.create_meshset(MESHSET_SET,out_meshset_tris); CHKERRG(rval);
  Range tris;
  rval = moab.get_adjacencies(prisms,2,false,tris,moab::Interface::UNION); CHKERRG(rval);
  std::cerr << tris.size() << " : " << prisms.size() << std::endl;
  rval = moab.add_entities(out_meshset_tris,tris); CHKERRG(rval);
  rval = moab.write_file("out_tris.vtk","VTK","",&out_meshset_tris,1); CHKERRG(rval);


  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  MoFEM::Core::Finalize();

}
