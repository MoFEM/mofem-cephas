/** \file dm_mofem.cpp

  \brief Atom test for Data Manager Interface in MoFEM

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

static const bool debug = false;

int main(int argc, char *argv[]) {




    //initialize petsc
    MoFEM::Core::Initialize(&argc,&argv,(char *)0,help);

    try {

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

    //register new dm type, i.e. mofem
    DMType dm_name = "DMMOFEM";
    ierr = DMRegister_MoFEM(dm_name); CHKERRG(ierr);
    DMType dm_name_sub = "DMSUB";
    ierr = DMRegister_MoFEM(dm_name_sub); CHKERRG(ierr);

    //craete dm instance
    DM dm;
    ierr = DMCreate(PETSC_COMM_WORLD,&dm);CHKERRG(ierr);
    ierr = DMSetType(dm,dm_name);CHKERRG(ierr);

    DM subdm0,subdm1;
    ierr = DMCreate(PETSC_COMM_WORLD,&subdm0);CHKERRG(ierr);
    ierr = DMSetType(subdm0,dm_name_sub);CHKERRG(ierr);
    ierr = DMCreate(PETSC_COMM_WORLD,&subdm1);CHKERRG(ierr);
    ierr = DMSetType(subdm1,dm_name_sub);CHKERRG(ierr);

    //read mesh and create moab and mofem datastrutures
    moab::Core mb_instance;
    moab::Interface& moab = mb_instance;
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
    const char *option;
    option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);

    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    EntityHandle root_set = moab.get_root_set();
    //add all entities to database, all of them will be used
    BitRefLevel bit_level0;
    bit_level0.set(0);
    ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(root_set,3,bit_level0); CHKERRG(ierr);
    //define & build field
    const int field_rank = 1;
    ierr = m_field.add_field("FIELD0",H1,AINSWORTH_LEGENDRE_BASE,field_rank); CHKERRG(ierr);
    //add entities to field
    ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"FIELD0"); CHKERRG(ierr);
    //set app. order
    ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD0",1); CHKERRG(ierr);
    ierr = m_field.add_field("FIELD1",H1,AINSWORTH_LEGENDRE_BASE,field_rank); CHKERRG(ierr);
    //add entities to field
    ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"FIELD1"); CHKERRG(ierr);
    //set app. order
    ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRG(ierr);

    //build data structures for fields
    ierr = m_field.build_fields(); CHKERRG(ierr);

    //define & build finite elements
    ierr = m_field.add_finite_element("FE00"); CHKERRG(ierr);
    //Define rows/cols and element data
    ierr = m_field.modify_finite_element_add_field_row("FE00","FIELD0"); CHKERRG(ierr);
    ierr = m_field.modify_finite_element_add_field_col("FE00","FIELD0"); CHKERRG(ierr);
    ierr = m_field.modify_finite_element_add_field_data("FE00","FIELD0"); CHKERRG(ierr);
    //add entities to finite element
    ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"FE00"); CHKERRG(ierr);

    //define & build finite elements
    ierr = m_field.add_finite_element("FE11"); CHKERRG(ierr);
    //Define rows/cols and element data
    ierr = m_field.modify_finite_element_add_field_row("FE11","FIELD1"); CHKERRG(ierr);
    ierr = m_field.modify_finite_element_add_field_col("FE11","FIELD1"); CHKERRG(ierr);
    ierr = m_field.modify_finite_element_add_field_data("FE11","FIELD1"); CHKERRG(ierr);
    //add entities to finite element
    ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"FE11"); CHKERRG(ierr);

    //define & build finite elements
    ierr = m_field.add_finite_element("FE01"); CHKERRG(ierr);
    //Define rows/cols and element data
    ierr = m_field.modify_finite_element_add_field_row("FE01","FIELD0"); CHKERRG(ierr);
    ierr = m_field.modify_finite_element_add_field_col("FE01","FIELD1"); CHKERRG(ierr);
    ierr = m_field.modify_finite_element_add_field_data("FE01","FIELD0"); CHKERRG(ierr);
    ierr = m_field.modify_finite_element_add_field_data("FE01","FIELD1"); CHKERRG(ierr);
    //add entities to finite element
    ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"FE01"); CHKERRG(ierr);

    //build finite elemnts
    ierr = m_field.build_finite_elements(); CHKERRG(ierr);
    //build adjacencies
    ierr = m_field.build_adjacencies(bit_level0); CHKERRG(ierr);

    //set dm data structure which created mofem data structures
    ierr = DMMoFEMCreateMoFEM(dm,&m_field,dm_name,bit_level0); CHKERRG(ierr);
    ierr = DMSetFromOptions(dm); CHKERRG(ierr);
    ierr = DMMoFEMAddElement(dm,"FE00"); CHKERRG(ierr);
    ierr = DMMoFEMAddElement(dm,"FE11"); CHKERRG(ierr);
    ierr = DMMoFEMAddElement(dm,"FE01"); CHKERRG(ierr);
    ierr = DMSetUp(dm); CHKERRG(ierr);
    ierr = m_field.partition_check_matrix_fill_in(dm_name,-1,-1,1); CHKERRG(ierr);

    int nf;
    char **field_names;
    IS *fields;
    ierr = DMCreateFieldIS(dm,&nf,&field_names,&fields); CHKERRG(ierr);
    for(int f = 0;f!=nf;f++) {
      PetscPrintf(PETSC_COMM_WORLD,"%d field %s\n",f,field_names[f]);
      ierr = PetscFree(field_names[f]); CHKERRG(ierr);
      ierr = ISDestroy(&(fields[f])); CHKERRG(ierr);
    }
    ierr = PetscFree(field_names); CHKERRG(ierr);
    ierr = PetscFree(fields); CHKERRG(ierr);

    ierr = DMMoFEMCreateSubDM(subdm0,dm,"SUB0"); CHKERRG(ierr);
    ierr = DMMoFEMSetSquareProblem(subdm0,PETSC_TRUE); CHKERRG(ierr);
    ierr = DMMoFEMAddElement(subdm0,"FE11"); CHKERRG(ierr);
    ierr = DMMoFEMAddSubFieldRow(subdm0,"FIELD1"); CHKERRG(ierr);
    ierr = DMMoFEMAddSubFieldCol(subdm0,"FIELD1"); CHKERRG(ierr);
    ierr = DMSetUp(subdm0); CHKERRG(ierr);
    ierr = m_field.partition_check_matrix_fill_in("SUB0",-1,-1,1); CHKERRG(ierr);
    if(debug) {
      Mat A;
      ierr = DMCreateMatrix(subdm0,&A); CHKERRG(ierr);
      MatView(A,PETSC_VIEWER_DRAW_WORLD);
      std::string wait;
      std::cin >> wait;
      ierr = MatDestroy(&A); CHKERRG(ierr);
    }

    ierr = DMMoFEMCreateSubDM(subdm1,dm,"SUB1"); CHKERRG(ierr);
    ierr = DMMoFEMSetSquareProblem(subdm1,PETSC_FALSE); CHKERRG(ierr);
    ierr = DMMoFEMAddElement(subdm1,"FE01"); CHKERRG(ierr);
    ierr = DMMoFEMAddSubFieldRow(subdm1,"FIELD0"); CHKERRG(ierr);
    ierr = DMMoFEMAddSubFieldCol(subdm1,"FIELD1"); CHKERRG(ierr);
    ierr = DMSetUp(subdm1); CHKERRG(ierr);
    ierr = m_field.partition_check_matrix_fill_in("SUB1",-1,-1,1); CHKERRG(ierr);

    if(debug) {
      Mat B;
      ierr = DMCreateMatrix(subdm1,&B); CHKERRG(ierr);
      MatView(B,PETSC_VIEWER_DRAW_WORLD);
      std::string wait;
      std::cin >> wait;
      ierr = MatDestroy(&B); CHKERRG(ierr);
    }

    //destry dm
    ierr = DMDestroy(&dm); CHKERRG(ierr);
    ierr = DMDestroy(&subdm0); CHKERRG(ierr);
    ierr = DMDestroy(&subdm1); CHKERRG(ierr);

    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    }

    //finish work cleaning memory, getting statistics, ect.
    ierr = MoFEM::Core::Finalize(); CHKERRG(ierr);

    return 0;

  }
