/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 *

 * This file is part of MoFEM.
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. 

 */

#include <tetgen.h>
#ifdef REAL
  #undef REAL
#endif

#include <MoFEM.hpp>
#include <TetGenInterface.hpp>
#include <BitLevelCoupler.hpp>

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

  //mapping between MoAB and TetGen
  tetgenio in,out;
  map<EntityHandle,unsigned long> moab_tetgen_map;
  map<unsigned long,EntityHandle> tetgen_moab_map;

  //get TetGen interface
  TetGenInterface *tetgen_iface;
  ierr = m_field.query_interface(tetgen_iface); CHKERRQ(ierr);

  //set MoAB nodes to TetGen data structure
  ierr = tetgen_iface->inData(nodes,in,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);
  map<int,Range> types_ents;
  types_ents[TetGenInterface::RIDGEVERTEX].merge(nodes);
  ierr = tetgen_iface->setGeomData(in,moab_tetgen_map,tetgen_moab_map,types_ents); CHKERRQ(ierr);

  char switches[] = ""; // TetGen switches, look to TegGen user manual
  ierr = tetgen_iface->tetRahedralize(switches,in,out); CHKERRQ(ierr); // make a TetGen mesh
  BitRefLevel bit_level1;
  bit_level1.set(1);
  //get mesh from TetGen and set bit level to 1
  ierr = tetgen_iface->outData(in,out,moab_tetgen_map,tetgen_moab_map,bit_level1,false,true); CHKERRQ(ierr);

  //save results
  EntityHandle meshset_level1;
  rval = moab.create_meshset(MESHSET_SET,meshset_level1); CHKERR_PETSC(rval);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level1,BitRefLevel().set(),MBTET,meshset_level1); CHKERRQ(ierr);
  if(debug) rval = moab.write_file("level1.vtk","VTK","",&meshset_level1,1); CHKERR_PETSC(rval);

  //clean data
  in.deinitialize();
  out.deinitialize();
  //initialize TetGen data structures
  in.initialize();
  out.initialize();

  //claer data stratus used by MoFEM
  moab_tetgen_map.clear();
  tetgen_moab_map.clear();

  //get mesh skin
  Skinner skin(&moab);
  Range outer_surface_skin;
  rval = skin.find_skin(0,tets,false,outer_surface_skin); CHKERR(rval);

  //set nodes to TetGen data structures
  ierr = tetgen_iface->inData(nodes,in,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);
  
  
  Range side_set_faces;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,sit)) {
    Range faces;
    rval = moab.get_entities_by_type(sit->meshset,MBTRI,faces,true); CHKERR_PETSC(rval);
    side_set_faces.merge(faces);
  }

  outer_surface_skin = subtract(outer_surface_skin,side_set_faces);
  vector<vector<Range> > sorted_outer_surface_skin;
  int nb_ents = outer_surface_skin.size();
  //sort surface elements in planar groups
  ierr = tetgen_iface->groupRegion_Triangle(outer_surface_skin,sorted_outer_surface_skin,1e-10); CHKERRQ(ierr);
  if(debug>0) {
    cout << " number of triangle entities " << nb_ents
      << " number of faces disjoint regions " << sorted_outer_surface_skin.size() << endl;
    for(unsigned int vv = 0;vv<sorted_outer_surface_skin.size();vv++) {
      cout << "\tnb of disjoint region " << vv 
	<< " nb of no-planar subregions " << sorted_outer_surface_skin[vv].size() << endl;
      for(unsigned int vvv = 0;vvv<sorted_outer_surface_skin[vv].size();vvv++) {
	cout << "\t\tnb. of subregion " << vvv << " nb. elements in subregion " << sorted_outer_surface_skin[vv][vvv].size() << endl;
      }
    }
  }

  vector<pair<Range,int> > markers; // set markers to surface elements
  vector<vector<Range> >::iterator vit = sorted_outer_surface_skin.begin();
  for(;vit!=sorted_outer_surface_skin.end();vit++) {
    vector<Range>::iterator viit = vit->begin();
    for(;viit!=vit->end();viit++) {
      Range polygons;
      ierr = tetgen_iface->makePolygonFacet(*viit,polygons); CHKERRQ(ierr);
      Range aa;
      rval = moab.get_connectivity(*viit,aa,true); CHKERR_PETSC(rval);
      Range viit_edges;
      rval = m_field.get_moab().get_adjacencies(
        *viit,1,false,viit_edges,Interface::UNION); CHKERR_PETSC(rval);
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,sit)) {
	Range faces;
	rval = moab.get_entities_by_type(sit->meshset,MBTRI,faces,true); CHKERR_PETSC(rval);
	Range faces_edges;
	rval = m_field.get_moab().get_adjacencies(
	  faces,1,false,faces_edges,Interface::UNION); CHKERR_PETSC(rval);
	aa.merge(intersect(faces_edges,viit_edges));
      }
      markers.push_back(pair<Range,int>(unite(polygons,aa),-1));
      //markers.push_back(pair<Range,int>(polygons,-1));
    }
  }
  
  // set markers to side Cubit sets
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,sit)) {
    int id = sit->get_msId();
    Range faces;
    rval = moab.get_entities_by_type(sit->meshset,MBTRI,faces,true); CHKERR_PETSC(rval);
    markers.push_back(pair<Range,int>(faces,id));
  }

  //ierr = tetgen_iface->inData(nodes,in,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);
  //set face markers
  ierr = tetgen_iface->setFaceData(markers,in,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);

  vector<pair<EntityHandle,int> > regions; // list of regions
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,bit)) {
    int id = bit->get_msId();
    Range tets;
    rval = moab.get_entities_by_type(bit->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
    regions.push_back(pair<EntityHandle,int>(*tets.begin(),-id));
  }
  //set volume regions
  ierr = tetgen_iface->setReginData(regions,in);  CHKERRQ(ierr);
  
  //print mesh in tetgeb format
  if(debug>0) {
    char tetgen_in_file_name[] = "in";
    in.save_nodes(tetgen_in_file_name);
    in.save_elements(tetgen_in_file_name);
    in.save_faces(tetgen_in_file_name);
    in.save_edges(tetgen_in_file_name);
    in.save_poly(tetgen_in_file_name);
  }

  //make TetGen mesh
  char switches2[] = "pYA"; // TetGen line command switches, look to TetGen manual for details
  ierr = tetgen_iface->tetRahedralize(switches2,in,out); CHKERRQ(ierr); // make a mesh
  BitRefLevel bit_level2;
  bit_level2.set(2);
  //get mesh form TetGen and store it on second bit refined level
  ierr = tetgen_iface->outData(in,out,moab_tetgen_map,tetgen_moab_map,bit_level2,false,true); CHKERRQ(ierr);
  //set markers to triangles
  ierr = tetgen_iface->getTriangleMarkers(tetgen_moab_map,out); CHKERRQ(ierr);
  //set refions to triangles
  ierr = tetgen_iface->getReginData(tetgen_moab_map,out); CHKERRQ(ierr);

  //char tetgen_out_file_name[] = "out";
  //out.save_elements(tetgen_out_file_name);

  //post-process results, save a mesh and skin
  EntityHandle meshset_level2;
  rval = moab.create_meshset(MESHSET_SET,meshset_level2); CHKERR_PETSC(rval);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level2,BitRefLevel().set(),MBTET,meshset_level2); CHKERRQ(ierr);
  if(debug) rval = moab.write_file("level2.vtk","VTK","",&meshset_level2,1); CHKERR_PETSC(rval);
  EntityHandle meshset_skin_level2;
  rval = moab.create_meshset(MESHSET_SET,meshset_skin_level2); CHKERR_PETSC(rval);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level2,BitRefLevel().set(),MBTRI,meshset_skin_level2); CHKERRQ(ierr);
  if(debug) rval = moab.write_file("level_skin2.vtk","VTK","",&meshset_skin_level2,1); CHKERR_PETSC(rval);

  //test BitLevelCoupler
  BitLevelCouplerInterface *bit_ref_copuler_ptr;
  ierr = m_field.query_interface(bit_ref_copuler_ptr); CHKERRQ(ierr);

  Range children;
  ierr = m_field.get_entities_by_ref_level(bit_level1,BitRefLevel().set(),children); CHKERRQ(ierr);
  if(children.empty()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"it should not be empty");
  }
  ierr = bit_ref_copuler_ptr->buidlAdjacenciesVerticesOnTets(bit_level0,children,true,1e-10,1e-6,true,2); CHKERRQ(ierr);
  ierr = bit_ref_copuler_ptr->buidlAdjacenciesEdgesFacesVolumes(bit_level0,children,true,2); CHKERRQ(ierr);

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  PetscFinalize();

}

