/** \file reading_med.cpp

  \brief Partition mesh and configuring blocksets

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

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  
  

  PetscInitialize(&argc,&argv,(char *)0,help);


  //global variables
  char mesh_file_name[255];
  PetscBool flg_file = PETSC_FALSE;
  PetscBool flg_n_part = PETSC_FALSE;
  PetscInt n_partas = 1;
  PetscBool creare_lower_dim_ents = PETSC_TRUE;

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","none","none"); CHKERRQ(ierr);
  ierr = PetscOptionsString(
    "-my_file",
    "mesh file name","",
    "mesh.h5m",mesh_file_name,255,&flg_file
  ); CHKERRQ(ierr);

  ierr = PetscOptionsInt(
    "-my_nparts",
    "number of parts","",
    1,&n_partas,&flg_n_part
  ); CHKERRQ(ierr);

  ierr = PetscOptionsBool(
    "-my_create_lower_dim_ents",
    "if tru create lower dimension entitities","",
    creare_lower_dim_ents,&creare_lower_dim_ents,NULL
  ); CHKERRQ(ierr);


  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);

  //Create MoFEM  database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  if(flg_file != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  if(flg_n_part != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR partitioning number not given");
  }

  MeshsetsManager *meshsets_interface_ptr;
  ierr = m_field.getInterface(meshsets_interface_ptr); CHKERRQ(ierr);
  ierr = meshsets_interface_ptr->setMeshsetFromFile(); CHKERRQ(ierr);

  for(
    CubitMeshSet_multiIndex::iterator cit = meshsets_interface_ptr->getBegin();
    cit!=meshsets_interface_ptr->getEnd(); cit++
  ) {
    std::cout << *cit << endl;
  }

  {
    Range ents3d;
    rval = moab.get_entities_by_dimension(0,3,ents3d,false); CHKERRQ_MOAB(rval);
    if(creare_lower_dim_ents) {
      Range faces;
      rval = moab.get_adjacencies(ents3d,2,true,faces,moab::Interface::UNION); CHKERRQ_MOAB(rval);
      Range edges;
      rval = moab.get_adjacencies(ents3d,1,true,edges,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    }
    ProblemsManager *prb_mng_ptr;
    ierr = m_field.getInterface(prb_mng_ptr); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionMesh(ents3d,3,2,n_partas); CHKERRQ(ierr);
  }

  rval = moab.write_file("out.h5m"); CHKERRQ_MOAB(rval);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}
