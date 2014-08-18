/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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


#include <tetgen.h>
#ifdef REAL
  #undef REAL
#endif

#include <MoFEM.hpp>
#include <TetGenInterface.hpp>

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //Read parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) databas
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
    
  Range nodes;
  rval = moab.get_entities_by_type(0,MBVERTEX,nodes,false); CHKERR_PETSC(rval);

  tetgenio in,out;
  map<EntityHandle,unsigned long> moab_tetgen_map;
  map<unsigned long,EntityHandle> tetgen_moab_map;

  TetGenInterface *tetgen_iface;
  ierr = m_field.query_interface(tetgen_iface); CHKERRQ(ierr);

  ierr = tetgen_iface->inData(nodes,in,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);
  char switches[] = "";
  ierr = tetgen_iface->tetRahedralize(switches,in,out); CHKERRQ(ierr);
  Range ents_out;
  ierr = tetgen_iface->outData(ents_out,in,out,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);

  BitRefLevel bit_level1;
  bit_level1.set(1);
  ierr = m_field.seed_ref_level_3D(ents_out,bit_level1); CHKERRQ(ierr);
  EntityHandle meshset_level1;
  rval = moab.create_meshset(MESHSET_SET,meshset_level1); CHKERR_PETSC(rval);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level1,BitRefLevel().set(),MBTET,meshset_level1); CHKERRQ(ierr);
  rval = moab.write_file("level1.vtk","VTK","",&meshset_level1,1); CHKERR_PETSC(rval);
        
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  PetscFinalize();

}

