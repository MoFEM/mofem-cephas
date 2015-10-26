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
#include <PrismsFromSurfaceInterface.hpp>

using namespace MoFEM;

MoABErrorCode rval;
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

    PrismsFromSurfaceInterface *prisms_from_surface_interface;
    ierr = m_field.query_interface(prisms_from_surface_interface); CHKERRQ(ierr);

    Range tris;
    rval = moab.get_entities_by_type(0,MBTRI,tris,false); CHKERR_PETSC(rval);
    Range prisms;
    ierr = prisms_from_surface_interface->createPrisms(tris,prisms); CHKERRQ(ierr);

    EntityHandle meshset;
    rval = moab.create_meshset(MESHSET_SET,meshset); CHKERR_PETSC(rval);
    rval = moab.add_entities(meshset,prisms); CHKERR_PETSC(rval);

    BitRefLevel bit_level0;
    bit_level0.set(0);
    ierr = m_field.seed_ref_level_3D(meshset,bit_level0); CHKERRQ(ierr);
    ierr = prisms_from_surface_interface->seedPrismsEntities(prisms,bit_level0); CHKERRQ(ierr);

    //Fields
    ierr = m_field.add_field("H1",H1,1); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_field_by_PRISMs(meshset,"H1",10); CHKERRQ(ierr);

    ierr = m_field.set_field_order(0,MBVERTEX,"H1",1,10); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBEDGE,"H1",2,10); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBTRI,"H1",3,10); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBQUAD,"H1",4,10); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBPRISM,"H1",6,10); CHKERRQ(ierr);
    ierr = m_field.build_fields(10); CHKERRQ(ierr);


    if(debug) {
      rval = moab.write_file("prism_mesh.vtk","VTK","",&meshset,1); CHKERR_PETSC(rval);
    }

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  PetscFinalize();
  return 0;

}
