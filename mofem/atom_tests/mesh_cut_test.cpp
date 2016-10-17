/** \file mesh_cut_test.cpp
  * \brief test for mesh cut interface
  *
  * \ingroup mesh_cut
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

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "testing mesh cut test\n\n";

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

  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    cout << *it << endl;
  }

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  // get cut mesh interface
  CutMeshInterface *cut_mesh;
  ierr = m_field.query_interface(cut_mesh); CHKERRQ(ierr);

  // get meshset manager interface
  MeshsetsManager *meshset_manager;
  ierr = m_field.query_interface(meshset_manager); CHKERRQ(ierr);

  // get surface entities
  Range surface;
  ierr = meshset_manager->getEntitiesByDimension(1,SIDESET,2,surface,true); CHKERRQ(ierr);
  cout << surface << endl;

  // get volume entities
  Range tets;
  ierr = moab.get_entities_by_dimension(0,3,tets,false); CHKERRQ(ierr);
  std::cout << tets << endl;

  // set mesh cutter entities
  ierr = cut_mesh->setSurface(surface); CHKERRQ(ierr);
  ierr = cut_mesh->setVolume(tets); CHKERRQ(ierr);

  // build tree
  ierr = cut_mesh->buildTree(); CHKERRQ(ierr);

  // find edges to cut
  ierr = cut_mesh->findToCut(0); CHKERRQ(ierr);

  // find tets at the front
  ierr = cut_mesh->findTetOnTheFront(5); CHKERRQ(ierr);

  // cut tets
  BitRefLevel bit_level1;
  bit_level1.set(1);
  #ifdef WITH_TETGEN
  ierr = cut_mesh->imprintFront(bit_level1,175,10); CHKERRQ(ierr);
  #endif

  // ierr = cut_mesh->cutTets(bit_level1); CHKERRQ(ierr);
  // ierr = cut_mesh->moveNodes(); CHKERRQ(ierr);

  EntityHandle meshset_vol;
  rval = moab.create_meshset(MESHSET_SET,meshset_vol); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(meshset_vol,tets); CHKERRQ_MOAB(rval);
  rval = moab.write_file("out_vol.vtk","VTK","",&meshset_vol,1); CHKERRQ_MOAB(rval);

  EntityHandle meshset_surf;
  rval = moab.create_meshset(MESHSET_SET,meshset_surf); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(meshset_surf,surface); CHKERRQ_MOAB(rval);
  rval = moab.write_file("out_surf.vtk","VTK","",&meshset_surf,1); CHKERRQ_MOAB(rval);

  EntityHandle meshset_cut_edges;
  rval = moab.create_meshset(MESHSET_SET,meshset_cut_edges); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(meshset_cut_edges,cut_mesh->getCutEdges()); CHKERRQ_MOAB(rval);
  rval = moab.write_file("out_cut_edges.vtk","VTK","",&meshset_cut_edges,1); CHKERRQ_MOAB(rval);

  EntityHandle meshset_cut_vol;
  rval = moab.create_meshset(MESHSET_SET,meshset_cut_vol); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(meshset_cut_vol,cut_mesh->getCutVolumes()); CHKERRQ_MOAB(rval);
  rval = moab.write_file("out_cut_vol.vtk","VTK","",&meshset_cut_vol,1); CHKERRQ_MOAB(rval);

  EntityHandle meshset_cut_tets;
  rval = moab.create_meshset(MESHSET_SET,meshset_cut_tets); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(meshset_cut_tets,cut_mesh->getNewCutVolumes()); CHKERRQ_MOAB(rval);
  rval = moab.write_file("out_cut_tets.vtk","VTK","",&meshset_cut_tets,1); CHKERRQ_MOAB(rval);

  EntityHandle meshset_new_faces;
  rval = moab.create_meshset(MESHSET_SET,meshset_new_faces); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(meshset_new_faces,cut_mesh->getNewCutSurfaces()); CHKERRQ_MOAB(rval);
  rval = moab.write_file("out_cut_new_faces.vtk","VTK","",&meshset_new_faces,1); CHKERRQ_MOAB(rval);

  EntityHandle meshset_front_tets;
  rval = moab.create_meshset(MESHSET_SET,meshset_front_tets); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(meshset_front_tets,cut_mesh->getFrontTets()); CHKERRQ_MOAB(rval);
  rval = moab.write_file("out_cut_front_tets.vtk","VTK","",&meshset_front_tets,1); CHKERRQ_MOAB(rval);


  PetscFinalize();

  return 0;
}
