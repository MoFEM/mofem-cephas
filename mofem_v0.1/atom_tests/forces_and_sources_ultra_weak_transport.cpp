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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"

#include "ForcesAndSurcesCore.hpp"
#include "TsCtx.hpp"
#include "UltaWeakTransportElement.hpp"
#include "DirichletBC.hpp"

#include "FEM.h"

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

#include <petscksp.h>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
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

  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& m_field = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("FLUXES",HDIV,1); CHKERRQ(ierr);
  ierr = m_field.add_field("VALUES",L2,1); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"FLUXES"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"VALUES"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 1;
  ierr = m_field.set_field_order(root_set,MBTET,"FLUXES",order+1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FLUXES",order+1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTET,"VALUES",order); CHKERRQ(ierr);

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  //finite elements

  /** thefine sources and other stuff 
    *
    * UltraWeakTransportElement is a class collecting functons, opertors and
    * data for ultra week implementation of transport element. See there to
    * learn how elements are created or how operators look like.
    *
    * Some methods in UltraWeakTransportElement are abstract, f.e. user need to
    * implement own surce therm.
    *
    */
  struct MyUltraWeakFE: public UltraWeakTransportElement {
   
    MyUltraWeakFE(FieldInterface &m_field): UltraWeakTransportElement(m_field) {};
 
    PetscErrorCode getFlux(const double x,const double y,const double z,double &flux) {
      PetscFunctionBegin;
      double d = sqrt(x*x+y*y+z*z);
      flux = 1-pow(d,5./4.);
      PetscFunctionReturn(0);
    }

  };

  MyUltraWeakFE ufe(m_field);
  ierr = ufe.addFiniteElements("FLUXES","VALUES"); CHKERRQ(ierr);

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("ULTRAWEAK"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("ULTRAWEAK",bit_level0); CHKERRQ(ierr);
  
  ierr = m_field.modify_problem_add_finite_element("ULTRAWEAK","ULTRAWEAK_FLUXFLUX"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ULTRAWEAK","ULTRAWEAK_FLUXVALUE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ULTRAWEAK","ULTRAWEAK_VALUEFLUX"); CHKERRQ(ierr);

  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //mesh partitioning 
  //partition
  ierr = m_field.simple_partition_problem("ULTRAWEAK"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("ULTRAWEAK"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("ULTRAWEAK"); CHKERRQ(ierr);

  Mat Aij;
  ierr = m_field.MatCreateMPIAIJWithArrays("ULTRAWEAK",&Aij); CHKERRQ(ierr);

  MatView(Aij,PETSC_VIEWER_DRAW_WORLD);
  std::string wait;
  std::cin >> wait;

  ierr = MatDestroy(&Aij); CHKERRQ(ierr);


  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


