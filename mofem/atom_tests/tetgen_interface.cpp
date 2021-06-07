/** \file tetgen_interface.cpp
  * \brief test for TetGen interface
  *
  * \ingroup mesh_tetgen
  */

/*
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

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";
static int debug = 1;

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc,&argv,(char *)0,help);

  try {

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //Read parameters from line command
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

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);

  //Create MoFEM (Joseph) databas
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(0,3,bit_level0); CHKERRG(ierr);

  Range nodes,tets;
  rval = moab.get_entities_by_type(0,MBVERTEX,nodes,false); CHKERRG(rval);
  rval = moab.get_entities_by_type(0,MBTET,tets,false); CHKERRG(rval);

  //mapping between MoAB and TetGen
  tetgenio in,out;
  std::map<EntityHandle,unsigned long> moab_tetgen_map;
  std::map<unsigned long,EntityHandle> tetgen_moab_map;

  //get TetGen interface
  TetGenInterface *tetgen_iface;
  ierr = m_field.getInterface(tetgen_iface); CHKERRG(ierr);

  //set MoAB nodes to TetGen data structure
  ierr = tetgen_iface->inData(nodes,in,moab_tetgen_map,tetgen_moab_map); CHKERRG(ierr);
  std::map<int,Range> types_ents;
  types_ents[TetGenInterface::RIDGEVERTEX].merge(nodes);
  ierr = tetgen_iface->setGeomData(in,moab_tetgen_map,tetgen_moab_map,types_ents); CHKERRG(ierr);

  char switches[] = ""; // TetGen switches, look to TegGen user manual
  ierr = tetgen_iface->tetRahedralize(switches,in,out); CHKERRG(ierr); // make a TetGen mesh
  BitRefLevel bit_level1;
  bit_level1.set(1);
  //get mesh from TetGen and set bit level to 1
  ierr = tetgen_iface->outData(in,out,moab_tetgen_map,tetgen_moab_map,bit_level1,false,true); CHKERRG(ierr);

  //save results
  EntityHandle meshset_level1;
  rval = moab.create_meshset(MESHSET_SET,meshset_level1); CHKERRG(rval);
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(bit_level1,BitRefLevel().set(),MBTET,meshset_level1); CHKERRG(ierr);
  if(debug) rval = moab.write_file("level1.vtk","VTK","",&meshset_level1,1); CHKERRG(rval);

  //clean data
  in.deinitialize();
  out.deinitialize();
  //initialize TetGen data structures
  in.initialize();
  out.initialize();

  //clear data stratus used by MoFEM
  moab_tetgen_map.clear();
  tetgen_moab_map.clear();

  //get mesh skin
  Skinner skin(&moab);
  Range outer_surface_skin;
  rval = skin.find_skin(0,tets,false,outer_surface_skin); CHKERRG(rval);

  //set nodes to TetGen data structures
  ierr = tetgen_iface->inData(nodes,in,moab_tetgen_map,tetgen_moab_map); CHKERRG(ierr);


  Range side_set_faces;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,sit)) {
    Range faces;
    rval = moab.get_entities_by_type(sit->meshset,MBTRI,faces,true); CHKERRG(rval);
    side_set_faces.merge(faces);
  }

  outer_surface_skin = subtract(outer_surface_skin,side_set_faces);
  std::vector<std::vector<Range> > sorted_outer_surface_skin;
  int nb_ents = outer_surface_skin.size();
  //sort surface elements in planar groups
  ierr = tetgen_iface->groupRegion_Triangle(outer_surface_skin,sorted_outer_surface_skin,1e-10); CHKERRG(ierr);
  if(debug>0) {
    std::cout << " number of triangle entities " << nb_ents
    << " number of faces disjoint regions " << sorted_outer_surface_skin.size() << std::endl;
    for(unsigned int vv = 0;vv<sorted_outer_surface_skin.size();vv++) {
      std::cout << "\tnb of disjoint region " << vv
      << " nb of no-planar subregions " << sorted_outer_surface_skin[vv].size() << std::endl;
      for(unsigned int vvv = 0;vvv<sorted_outer_surface_skin[vv].size();vvv++) {
        std::cout << "\t\tnb. of subregion " << vvv << " nb. elements in subregion " << sorted_outer_surface_skin[vv][vvv].size() << std::endl;
      }
    }
  }

  std::vector<std::pair<Range,int> > markers; // set markers to surface elements
  std::vector<std::vector<Range> >::iterator vit = sorted_outer_surface_skin.begin();
  for(;vit!=sorted_outer_surface_skin.end();vit++) {
    std::vector<Range>::iterator viit = vit->begin();
    for(;viit!=vit->end();viit++) {
      Range polygons;
      ierr = tetgen_iface->makePolygonFacet(*viit,polygons); CHKERRG(ierr);
      Range aa;
      rval = moab.get_connectivity(*viit,aa,true); CHKERRG(rval);
      Range viit_edges;
      rval = m_field.get_moab().get_adjacencies(
        *viit,1,false,viit_edges,moab::Interface::UNION); CHKERRG(rval);
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,sit)) {
	      Range faces;
	      rval = moab.get_entities_by_type(sit->meshset,MBTRI,faces,true); CHKERRG(rval);
	      Range faces_edges;
	      rval = m_field.get_moab().get_adjacencies(
	      faces,1,false,faces_edges,moab::Interface::UNION); CHKERRG(rval);
	      aa.merge(intersect(faces_edges,viit_edges));
      }
      markers.push_back(std::pair<Range,int>(unite(polygons,aa),-1));
      //markers.push_back(std::pair<Range,int>(polygons,-1));
    }
  }

  // set markers to side Cubit sets
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,sit)) {
    int id = sit->getMeshsetId();
    Range faces;
    rval = moab.get_entities_by_type(sit->meshset,MBTRI,faces,true); CHKERRG(rval);
    markers.push_back(std::pair<Range,int>(faces,id));
  }

  //ierr = tetgen_iface->inData(nodes,in,moab_tetgen_map,tetgen_moab_map); CHKERRG(ierr);
  //set face markers
  ierr = tetgen_iface->setFaceData(markers,in,moab_tetgen_map,tetgen_moab_map); CHKERRG(ierr);

  std::vector<std::pair<EntityHandle,int> > regions; // list of regions
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,bit)) {
    int id = bit->getMeshsetId();
    Range tets;
    rval = moab.get_entities_by_type(bit->meshset,MBTET,tets,true); CHKERRG(rval);
    regions.push_back(std::pair<EntityHandle,int>(*tets.begin(),-id));
  }
  //set volume regions
  ierr = tetgen_iface->setRegionData(regions,in);  CHKERRG(ierr);

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
  ierr = tetgen_iface->tetRahedralize(switches2,in,out); CHKERRG(ierr); // make a mesh
  BitRefLevel bit_level2;
  bit_level2.set(2);
  //get mesh form TetGen and store it on second bit refined level
  ierr = tetgen_iface->outData(in,out,moab_tetgen_map,tetgen_moab_map,bit_level2,false,true); CHKERRG(ierr);
  //set markers to triangles
  ierr = tetgen_iface->getTriangleMarkers(tetgen_moab_map,out); CHKERRG(ierr);
  //set refions to triangles
  ierr = tetgen_iface->getRegionData(tetgen_moab_map,out); CHKERRG(ierr);

  //char tetgen_out_file_name[] = "out";
  //out.save_elements(tetgen_out_file_name);

  //post-process results, save a mesh and skin
  EntityHandle meshset_level2;
  rval = moab.create_meshset(MESHSET_SET,meshset_level2); CHKERRG(rval);
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(bit_level2,BitRefLevel().set(),MBTET,meshset_level2); CHKERRG(ierr);
  if(debug) rval = moab.write_file("level2.vtk","VTK","",&meshset_level2,1); CHKERRG(rval);
  EntityHandle meshset_skin_level2;
  rval = moab.create_meshset(MESHSET_SET,meshset_skin_level2); CHKERRG(rval);
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(bit_level2,BitRefLevel().set(),MBTRI,meshset_skin_level2); CHKERRG(ierr);
  if(debug) rval = moab.write_file("level_skin2.vtk","VTK","",&meshset_skin_level2,1); CHKERRG(rval);

  //test BitLevelCoupler
  BitLevelCoupler *bit_ref_coupler_ptr;
  ierr = m_field.getInterface(bit_ref_coupler_ptr); CHKERRG(ierr);

  Range children;
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(bit_level1,BitRefLevel().set(),children); CHKERRG(ierr);
  if(children.empty()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"it should not be empty");
  }
  ierr = bit_ref_coupler_ptr->buildAdjacenciesVerticesOnTets(bit_level0,children,true,1e-10,1e-6,true,2); CHKERRG(ierr);
  ierr = bit_ref_coupler_ptr->buildAdjacenciesEdgesFacesVolumes(bit_level0,children,true,2); CHKERRG(ierr);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  MoFEM::Core::Finalize();

}
