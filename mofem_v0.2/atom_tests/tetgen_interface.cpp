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

#include <moab/Skinner.hpp>

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
  
  Range nodes,tets;
  rval = moab.get_entities_by_type(0,MBVERTEX,nodes,false); CHKERR_PETSC(rval);
  rval = moab.get_entities_by_type(0,MBTET,tets,false); CHKERR_PETSC(rval);

  tetgenio in,out;
  map<EntityHandle,unsigned long> moab_tetgen_map;
  map<unsigned long,EntityHandle> tetgen_moab_map;

  TetGenInterface *tetgen_iface;
  ierr = m_field.query_interface(tetgen_iface); CHKERRQ(ierr);

  ierr = tetgen_iface->inData(nodes,in,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);
  char switches[] = "";
  ierr = tetgen_iface->tetRahedralize(switches,in,out); CHKERRQ(ierr);
  BitRefLevel bit_level1;
  bit_level1.set(1);
  ierr = tetgen_iface->outData(bit_level1,in,out,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);

  EntityHandle meshset_level1;
  rval = moab.create_meshset(MESHSET_SET,meshset_level1); CHKERR_PETSC(rval);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level1,BitRefLevel().set(),MBTET,meshset_level1); CHKERRQ(ierr);
  rval = moab.write_file("level1.vtk","VTK","",&meshset_level1,1); CHKERR_PETSC(rval);

  Skinner skin(&moab);
  Range outer_surface_skin;
  rval = skin.find_skin(0,tets,false,outer_surface_skin); CHKERR(rval);
  Range outer_surface_skin_nodes;
  moab.get_connectivity(outer_surface_skin,outer_surface_skin_nodes);
  outer_surface_skin.merge(outer_surface_skin_nodes);

  in.deinitialize();
  out.deinitialize();
  in.initialize();
  out.initialize();

  moab_tetgen_map.clear();
  tetgen_moab_map.clear();

  ierr = tetgen_iface->inData(outer_surface_skin,in,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);
  char switches2[] = "pY";
  ierr = tetgen_iface->tetRahedralize(switches2,in,out); CHKERRQ(ierr);
  BitRefLevel bit_level2;
  bit_level2.set(2);
  ierr = tetgen_iface->outData(bit_level2,in,out,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);

  EntityHandle meshset_level2;
  rval = moab.create_meshset(MESHSET_SET,meshset_level2); CHKERR_PETSC(rval);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level2,BitRefLevel().set(),MBTET,meshset_level2); CHKERRQ(ierr);
  rval = moab.write_file("level2.vtk","VTK","",&meshset_level2,1); CHKERR_PETSC(rval);

  in.deinitialize();
  out.deinitialize();
  in.initialize();
  out.initialize();

  moab_tetgen_map.clear();
  tetgen_moab_map.clear();

  tets.merge(nodes);
  tets.merge(outer_surface_skin);

  ierr = tetgen_iface->inData(tets,in,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);
  ierr = tetgen_iface->setFaceCubitSideSetMarkers(moab_tetgen_map,in); CHKERRQ(ierr);

  char switches3[] = "prq";
  ierr = tetgen_iface->tetRahedralize(switches3,in,out); CHKERRQ(ierr);
  BitRefLevel bit_level3;
  bit_level3.set(3);
  ierr = tetgen_iface->outData(bit_level3,in,out,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);
  ierr = tetgen_iface->getFaceCubitSideSetMarkers(tetgen_moab_map,out); CHKERRQ(ierr);

  Tag th;
  int def_marker = 0;
  rval = moab.tag_get_handle("MSIDSIDESET",1,MB_TYPE_INTEGER,th,MB_TAG_CREAT|MB_TAG_SPARSE,&def_marker); CHKERR_THROW(rval); 

  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,sit)) {

    int id = sit->get_msId();

    Range faces;
    rval = moab.get_entities_by_type(sit->meshset,MBTRI,faces,true); CHKERR_PETSC(rval);
    Range::iterator it = faces.begin();
    for(;it!=faces.end();it++) {
      rval = moab.tag_set_data(th,&*it,1,&id); CHKERR_PETSC(rval);
    }

  }

  EntityHandle meshset_level3;
  rval = moab.create_meshset(MESHSET_SET,meshset_level3); CHKERR_PETSC(rval);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level3,BitRefLevel().set(),MBTRI,meshset_level3); CHKERRQ(ierr);
  rval = moab.write_file("level3.vtk","VTK","",&meshset_level3,1); CHKERR_PETSC(rval);
        
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  PetscFinalize();

}

