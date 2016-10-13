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

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "testing mesh refinement algorithm\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_MOAB(rval);

  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;


  MeshRefinement *refine;
  ierr = m_field.query_interface(refine); CHKERRQ(ierr);

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  BitRefLevel bit_level1;
  bit_level1.set(1);

  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  // random mesh refinement
  EntityHandle meshset_ref_edges;
  rval = moab.create_meshset(MESHSET_SET,meshset_ref_edges); CHKERRQ_MOAB(rval);
  Range edges_to_refine;
  rval = moab.get_entities_by_type(meshset_level0,MBEDGE,edges_to_refine);  CHKERRQ_MOAB(rval);
  int ii = 0;
  for(Range::iterator eit = edges_to_refine.begin();
    eit!=edges_to_refine.end();eit++,ii++) {
    int numb = ii % 2;
    if(numb == 0) {
      ierr = moab.add_entities(meshset_ref_edges,&*eit,1); CHKERRQ(ierr);
    }
  }
  ierr = refine->add_verices_in_the_middel_of_edges(meshset_ref_edges,bit_level1); CHKERRQ(ierr);
  ierr = refine->refine_TET(meshset_level0,bit_level1); CHKERRQ(ierr);
  //PetscAttachDebugger ();
  //ierr = m_field.shift_right_bit_ref(1); CHKERRQ(ierr);

  std::ofstream myfile;
  myfile.open("mesh_refine.txt");

  EntityHandle out_meshset_tet;
  rval = moab.create_meshset(MESHSET_SET,out_meshset_tet); CHKERRQ_MOAB(rval);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level1,BitRefLevel().set(),MBTET,out_meshset_tet); CHKERRQ(ierr);
  Range tets;
  rval = moab.get_entities_by_handle(out_meshset_tet,tets); CHKERRQ_MOAB(rval);
  {
    int ii = 0;
    for(Range::iterator tit = tets.begin();tit!=tets.end();tit++) {
      int num_nodes;
      const EntityHandle* conn;
      rval = moab.get_connectivity(*tit,conn,num_nodes,true); CHKERRQ_MOAB(rval);

      for(int nn = 0;nn<num_nodes;nn++) {
        //cout << conn[nn] << " ";
        myfile << conn[nn] << " ";
      }
      //cout << std::endl;
      myfile << std::endl;
      if(ii>25) break;
    }
  }

  myfile.close();

  //rval = moab.write_file("out_mesh_refine.vtk","VTK","",&out_meshset_tet,1); CHKERRQ_MOAB(rval);

  BitLevelCouplerInterface *bit_ref_copuler_ptr;
  ierr = m_field.query_interface(bit_ref_copuler_ptr); CHKERRQ(ierr);

  Range children;
  ierr = m_field.get_entities_by_ref_level(bit_level1,BitRefLevel().set(),children); CHKERRQ(ierr);
  if(children.empty()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"it should not be empty");
  }
  bit_ref_copuler_ptr->vErify = true;
  ierr = bit_ref_copuler_ptr->buidlAdjacenciesEdgesFacesVolumes(bit_level0,children,true,2); CHKERRQ(ierr);

  //reset entities
  bit_ref_copuler_ptr->vErify = false;
  Range children_new;
  ierr = m_field.get_entities_by_ref_level(bit_level1,bit_level1,children_new); CHKERRQ(ierr);
  ierr = bit_ref_copuler_ptr->resetParents(children_new,true); CHKERRQ(ierr);

  ierr = bit_ref_copuler_ptr->buidlAdjacenciesVerticesOnTets(bit_level0,children,true,1e-10,1e-6,true,0); CHKERRQ(ierr);
  ierr = bit_ref_copuler_ptr->buidlAdjacenciesEdgesFacesVolumes(bit_level0,children,true,2); CHKERRQ(ierr);

  PetscFinalize();

  return 0;
}
