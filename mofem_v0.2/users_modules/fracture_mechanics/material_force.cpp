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

#include <MoFEM.hpp>
using namespace MoFEM;

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>
#include <ArcLengthTools.hpp>
#include <MatShellConstrainsByMarkAinsworth.hpp>

extern "C" {
  #include <complex_for_lazy.h>
}

#include <FEMethod_ComplexForLazy.hpp>
#include <FEMethod_DriverComplexForLazy.hpp>

#include <moab/Skinner.hpp>
#include <moab/AdaptiveKDTree.hpp>

#include <petsctime.h>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>
#include <PostProcNonLinearElasticityStresseOnRefindedMesh.hpp>

#include <FaceSplittingTool.hpp>
#include <ConfigurationalFractureMechanics.hpp>

using namespace ObosleteUsersModules;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  try {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm) 
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  BARRIER_RANK_END(pcomm) 

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  ConfigurationalFractureMechanics conf_prob(m_field);
  ierr = main_material_forces(m_field,conf_prob); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    rval = moab.write_file("out_material.h5m"); CHKERR_PETSC(rval);
  }

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = m_field.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = m_field.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  if(pcomm->rank()==0) {

    Range level_tris;
    ierr = m_field.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);

    Range SurfacesFaces;
    ierr = m_field.get_Cubit_msId_entities_by_dimension(102,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
    SurfacesFaces = intersect(SurfacesFaces,level_tris);
    Range CrackSurfacesFaces;
    ierr = m_field.get_Cubit_msId_entities_by_dimension(200,SIDESET,2,CrackSurfacesFaces,true); CHKERRQ(ierr);
    CrackSurfacesFaces = intersect(CrackSurfacesFaces,level_tris);
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
      int msId = it->get_msId();
      if((msId < 10200)||(msId >= 10300)) continue;
      Range SurfacesFaces_msId;
      ierr = m_field.get_Cubit_msId_entities_by_dimension(msId,SIDESET,2,SurfacesFaces_msId,true); CHKERRQ(ierr);
      SurfacesFaces_msId = intersect(SurfacesFaces_msId,level_tris);
      SurfacesFaces.insert(SurfacesFaces_msId.begin(),SurfacesFaces_msId.end());
    }

    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    rval = moab.add_entities(out_meshset,CrackSurfacesFaces); CHKERR_PETSC(rval);
    rval = moab.write_file("out_crack_surface.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.add_entities(out_meshset,SurfacesFaces); CHKERR_PETSC(rval);
    rval = moab.write_file("out_surface.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  return 0;
}
