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
namespace nglib {
#include <nglib.h>
}
using namespace nglib;
#include <NetGenInterface.hpp>

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";
static int debug = 1;

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //Read parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) databas
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  Skinner skin(&moab);
  // Define pointer to a new Netgen Mesh
  Ng_Mesh *mesh;
  // Define pointer to STL Geometry
  Ng_STL_Geometry *stl_geom;
  // Result of Netgen Operations
  Ng_Result ng_res;

  // Set the Meshing Parameters to be used
  Ng_Meshing_Parameters mp;
  mp.maxh = 0.1;
  mp.fineness = 0.5;
  mp.second_order = 0;

  NetGenInterface *netgen_iface;
  ierr = m_field.query_interface(netgen_iface); CHKERRQ(ierr);
  std::vector<EntityHandle> pts;
  EntityHandle meshset_out;

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  Range tets;
  rval = moab.get_entities_by_type(0,MBTET,tets,false); CHKERRQ_MOAB(rval);
  Range outer_surface_skin;
  rval = skin.find_skin(0,tets,false,outer_surface_skin); CHKERR_MOAB(rval);

  // Initialise the Netgen Core library
  Ng_Init();

  // Actually create the mesh structure
  mesh = Ng_NewMesh();
  // careate new geometry
  stl_geom = Ng_STL_NewGeometry();

  outer_surface_skin = unite(outer_surface_skin,tets);
  ierr = netgen_iface->stlSetSurfaceTriangles(stl_geom,outer_surface_skin,NULL,1);  CHKERRQ(ierr);

  cout << "Initialise the STL Geometry structure...." << std::endl;
  ng_res = Ng_STL_InitSTLGeometry(stl_geom);
  if(ng_res != NG_OK) {
    cout << "Error Initialising the STL Geometry....Aborting!!" << std::endl;
    return 1;
  }

  cout << "Start Edge Meshing...." << std::endl;
  ng_res = Ng_STL_MakeEdges(stl_geom, mesh, &mp);
  if(ng_res != NG_OK) {
    cout << "Error in Edge Meshing....Aborting!!" << std::endl;
    return 1;
  }


  cout << "Start Surface Meshing...." << std::endl;
  ng_res = Ng_STL_GenerateSurfaceMesh(stl_geom, mesh, &mp);
  if(ng_res != NG_OK) {
    cout << "Error in Surface Meshing....Aborting!!" << std::endl;
    return 1;
  }

  // volume mesh output
  int np;
  np = Ng_GetNP(mesh);
  cout << "Points: " << np << std::endl;

  int ne;
  ne = Ng_GetNSE(mesh);
  cout << "Elements Surface: " << ne << " (" << outer_surface_skin.size() << ") " << std::endl;

  ierr = netgen_iface->getPoints(mesh,pts); CHKERRQ(ierr);
  std::vector<EntityHandle> surface_elems;
  ierr = netgen_iface->getSurfaceElements(mesh,pts,surface_elems); CHKERRQ(ierr);

  rval = moab.create_meshset(MESHSET_SET,meshset_out); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(meshset_out,&*surface_elems.begin(),surface_elems.size()); CHKERR_MOAB(rval);
  if(debug) rval = moab.write_file("netgen1.vtk","VTK","",&meshset_out,1); CHKERRQ_MOAB(rval);
  rval = moab.delete_entities(&meshset_out,1); CHKERRQ_MOAB(rval);

  delete stl_geom;
  Ng_DeleteMesh(mesh);
  Ng_Exit();

  // Initialise the Netgen Core library
  Ng_Init();

  // Actually create the mesh structure
  mesh = Ng_NewMesh();

  ierr = netgen_iface->setPoints(mesh,pts); CHKERRQ(ierr);
  ierr = netgen_iface->setSurfaceElements(mesh,pts,surface_elems); CHKERRQ(ierr);

  np = Ng_GetNP(mesh);
  cout << "Points: " << np << std::endl;
  ne = Ng_GetNSE(mesh);
  cout << "Elements Surface: " << ne << " (" << outer_surface_skin.size() << ") " << std::endl;

  cout << "Start Volume Meshing...." << std::endl;
  ng_res = Ng_GenerateVolumeMesh (mesh, &mp);
  if(ng_res != NG_OK) {
    cout << "Error in Volume Meshing....Aborting!!" << std::endl;
    return 1;
  }
  cout << "Meshing successfully completed....!!" << std::endl;

  np = Ng_GetNP(mesh);
  cout << "Points: " << np << std::endl;
  ne = Ng_GetNE(mesh);
  cout << "Elements Volume: " << ne << std::endl;

  ierr = netgen_iface->getPoints(mesh,pts); CHKERRQ(ierr);
  std::vector<EntityHandle> volume_elems;
  ierr = netgen_iface->getVolumeElements(mesh,pts,volume_elems); CHKERRQ(ierr);

  rval = moab.create_meshset(MESHSET_SET,meshset_out); CHKERRQ_MOAB(rval);
  rval = moab.add_entities(meshset_out,&*volume_elems.begin(),volume_elems.size()); CHKERR_MOAB(rval);
  if(debug) rval = moab.write_file("netgen2.vtk","VTK","",&meshset_out,1); CHKERRQ_MOAB(rval);
  rval = moab.delete_entities(&meshset_out,1); CHKERRQ_MOAB(rval);

  Ng_DeleteMesh(mesh);
  Ng_Exit();

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFinalize();

}
