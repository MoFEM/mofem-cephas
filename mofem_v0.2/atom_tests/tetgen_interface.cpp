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
static int debug = 1;

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
  ierr = tetgen_iface->outData(in,out,moab_tetgen_map,tetgen_moab_map,bit_level1); CHKERRQ(ierr);

  EntityHandle meshset_level1;
  rval = moab.create_meshset(MESHSET_SET,meshset_level1); CHKERR_PETSC(rval);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level1,BitRefLevel().set(),MBTET,meshset_level1); CHKERRQ(ierr);
  if(debug) rval = moab.write_file("level1.vtk","VTK","",&meshset_level1,1); CHKERR_PETSC(rval);

  in.deinitialize();
  out.deinitialize();
  in.initialize();
  out.initialize();

  moab_tetgen_map.clear();
  tetgen_moab_map.clear();

  Skinner skin(&moab);
  Range outer_surface_skin;
  rval = skin.find_skin(0,tets,false,outer_surface_skin); CHKERR(rval);

  //ierr = tetgen_iface->inData(nodes,in,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);
  
  Range side_set_faces;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,sit)) {
    Range faces;
    rval = moab.get_entities_by_type(sit->meshset,MBTRI,faces,true); CHKERR_PETSC(rval);
    side_set_faces.merge(faces);
  }

  //Range surface_nodes;
  //rval = moab.get_connectivity(unite(side_set_faces,outer_surface_skin),surface_nodes,true); CHKERR_PETSC(rval);
  ierr = tetgen_iface->inData(nodes,in,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);

  vector<pair<Range,int> > markers;
  outer_surface_skin = subtract(outer_surface_skin,side_set_faces);
  Range::iterator it = outer_surface_skin.begin();
  for(;it!=outer_surface_skin.end();it++) {
    Range ent;
    ent.insert(*it);
    markers.push_back(pair<Range,int>(ent,-1));
  }
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,sit)) {
    int id = sit->get_msId();
    Range faces;
    rval = moab.get_entities_by_type(sit->meshset,MBTRI,faces,true); CHKERR_PETSC(rval);
    markers.push_back(pair<Range,int>(faces,id));
  }

  ierr = tetgen_iface->setFaceData(markers,in,moab_tetgen_map,tetgen_moab_map);  CHKERRQ(ierr);

  vector<pair<Range,int> > regions;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,bit)) {
    int id = bit->get_msId();
    Range tets;
    rval = moab.get_entities_by_type(bit->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
    Range one;
    one.insert(*tets.begin());
    regions.push_back(pair<Range,int>(one,id));
  }
  ierr = tetgen_iface->setReginData(regions,in);  CHKERRQ(ierr);
  
  //in.load_poly("bar2");
  //in.save_nodes("in");
  //in.save_poly("in");

  char switches2[] = "pYAz";
  ierr = tetgen_iface->tetRahedralize(switches2,in,out); CHKERRQ(ierr);
  BitRefLevel bit_level2;
  bit_level2.set(2);
  ierr = tetgen_iface->outData(in,out,moab_tetgen_map,tetgen_moab_map,bit_level2); CHKERRQ(ierr);
  ierr = tetgen_iface->getTiangleAttributes(tetgen_moab_map,out); CHKERRQ(ierr);
  ierr = tetgen_iface->getReginData(tetgen_moab_map,out); CHKERRQ(ierr);

  char tetgen_out_file_name[] = "out";
  out.save_elements(tetgen_out_file_name);

  EntityHandle meshset_level2;
  rval = moab.create_meshset(MESHSET_SET,meshset_level2); CHKERR_PETSC(rval);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level2,BitRefLevel().set(),MBTET,meshset_level2); CHKERRQ(ierr);
  if(debug) rval = moab.write_file("level2.vtk","VTK","",&meshset_level2,1); CHKERR_PETSC(rval);
  EntityHandle meshset_skin_level2;
  rval = moab.create_meshset(MESHSET_SET,meshset_skin_level2); CHKERR_PETSC(rval);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level2,BitRefLevel().set(),MBTRI,meshset_skin_level2); CHKERRQ(ierr);
  if(debug) rval = moab.write_file("level_skin2.vtk","VTK","",&meshset_skin_level2,1); CHKERR_PETSC(rval);

        
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  PetscFinalize();

}

