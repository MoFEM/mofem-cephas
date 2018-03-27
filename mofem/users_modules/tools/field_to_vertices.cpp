/** \file field_to_vertices.cpp
  \brief Field to vertices
  \example field_to_vertices.cpp

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

int main(int argc, char *argv[]) {
  PetscInitialize(&argc, &argv, (char *)0, help);

  try {

    // global variables
    char mesh_file_name[255];
    PetscBool flg_file = PETSC_FALSE;
    char field_name_param[255] = "RHO";
    CHKERR PetscOptionsBegin(PETSC_COMM_WORLD, "", "Field to vertices options", "none");
    CHKERR PetscOptionsString("-my_file", "mesh file name", "", "mesh.h5m",
                              mesh_file_name, 255, &flg_file);
    CHKERR PetscOptionsString("-my_field", "field name", "", "FIELD",
                              field_name_param, 255, PETSC_NULL);
    ierr = PetscOptionsEnd(); CHKERRG(ierr);

    std::string field_name(field_name_param);

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
    const char *option;
    option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // Create MoFEM  database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    if (flg_file != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
              "*** ERROR -my_file (MESH FILE NEEDED)");
    }
    // if (flg_list_of_sidesets != PETSC_TRUE) {
    //   SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
    //           "List of sidesets not given -my_side_sets ...");
    // }
    BitRefLevel bit_level0;
    bit_level0.set(0);
    ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

    ierr = m_field.build_fields();

    //ierr = m_field.list_fields(); CHKERRQ(ierr);

    // Projection10NodeCoordsOnField ents_method(m_field,"MESH_NODE_POSITIONS");
    // CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS",ents_method);

    if(m_field.check_field("rho")) {
      field_name = "rho";
    }

    bool field_flg = false;
    const Field_multiIndex *fields_ptr;
    CHKERR m_field.get_fields(&fields_ptr);
    for(auto field : (*fields_ptr)) {
      // cerr << field->getName() << "\n";
      // cerr << field->getNbOfCoeffs() << "\n";
      bool check_space = field->getSpace() == H1;
      if(field->getName() == field_name && check_space) field_flg = true;
    }
     if (!field_flg) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
              "*** ERROR -my_field (FIELD (in H1 space) is NOT FOUND)");
    }
    // const FieldEntity_multiIndex *field_ents;
    // CHKERR m_field.get_field_ents(field_ents);
    // auto low_it = field_ents.get<FieldName_mi_tag>.lower_bound("RHO");
    // auto hi_it = field_ents.get<FieldName_mi_tag>.upper_bound("RHO");
    // for(;low_it!=hi_it;++low_it) {
    //   if(low_it->get()->getEntType()==MBVERTEX) {
    //     double data = low_it->get()->getEntFieldData()[0];
    //     EntityHandle ent = low_it->get()->getEnt();
    //     CHKERR moab.tag_set_data(th,&ent,1,&data);

    //   }
    // }

    // PostProcVertexMethod ent_method(m_field, field_name.c_str());
    SaveVertexDofOnTag ent_method(m_field, field_name.c_str());

    CHKERR m_field.loop_dofs(field_name.c_str(),ent_method);
    // CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method);
    PetscPrintf(PETSC_COMM_WORLD, "\nDone. Saving files... \n");
            // Range nodes_tets;
            // CHKERR moab.get_entities_by_type(0, MBVERTEX, nodes_tets, true);
            // for ( auto it : nodes_tets){
            //   double var;
            //   moab.tag_get_data(ent_method.tH,&it,1,&var);
            //   cerr << "DATA on " << it << " is: " << var << '\n';
            // }

    //TODO: Higher order field mapping
    CHKERR m_field.getInterface<BitRefManager>()->writeBitLevelByType(bit_level0, BitRefLevel().set(), MBTET,
                                         "out_mesh.vtk", "VTK",
                                         "");
    CHKERR moab.write_file("out.h5m");
  }
  CATCH_ERRORS;

  CHKERR PetscFinalize();

  return 0;
}
