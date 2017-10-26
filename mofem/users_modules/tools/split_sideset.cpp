/** \file split_sideset.cpp
  \brief Split sidesets
  \example split_sideset.cpp

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
    PetscBool squash_bit_levels = PETSC_TRUE;
    PetscBool flg_list_of_sidesets = PETSC_FALSE;
    PetscBool output_vtk = PETSC_TRUE;

    int nb_sidesets = 10;
    int sidesets[nb_sidesets];

    ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "", "Split sides options", "none");
    CHKERRQ(ierr);
    ierr = PetscOptionsString("-my_file", "mesh file name", "", "mesh.h5m",
                              mesh_file_name, 255, &flg_file);
    CHKERRQ(ierr);

    ierr = PetscOptionsBool("-squash_bit_levels", "squahs bit levels", "",
                            squash_bit_levels, &squash_bit_levels, NULL);
    CHKERRQ(ierr);

    ierr = PetscOptionsIntArray("-side_sets", "get list of sidesets", "",
                                sidesets, &nb_sidesets, &flg_list_of_sidesets);
    CHKERRQ(ierr);
    ierr = PetscOptionsBool("-output_vtk", "if true outout vtk file", "",
                            output_vtk, &output_vtk, PETSC_NULL);
    CHKERRQ(ierr);
 

    ierr = PetscOptionsEnd();
    CHKERRQ(ierr);

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option);
    CHKERRQ_MOAB(rval);

    // Create MoFEM  database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    if (flg_file != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
              "*** ERROR -my_file (MESH FILE NEEDED)");
    }
    if (flg_list_of_sidesets != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
              "List of sidesets not given -my_side_sets ...");
    }

    // Get interface to meshsets manager
    MeshsetsManager *m_mng;
    ierr = m_field.getInterface(m_mng);
    CHKERRQ(ierr);
    // Get interface for splitting manager
    PrismInterface *interface_ptr;
    ierr = m_field.getInterface(interface_ptr);
    CHKERRQ(ierr);
    BitRefManager *bit_mng;
    ierr = m_field.getInterface(bit_mng);
    CHKERRQ(ierr);

    // Seed mesh with bit levels
    ierr = bit_mng->setBitRefLevelByDim(0, 3, BitRefLevel().set(0));
    CHKERRQ(ierr);
    std::vector<BitRefLevel> bit_levels;
    bit_levels.push_back(BitRefLevel().set(0));

    typedef CubitMeshSet_multiIndex::index<
      CubitMeshSets_mask_meshset_mi_tag>::type CMeshsetByType;
    typedef CMeshsetByType::iterator CMIteratorByType;
    typedef CubitMeshSet_multiIndex::index<
        Composite_Cubit_msId_And_MeshSetType_mi_tag>::type CMeshsetByIdType;
    typedef CMeshsetByIdType::iterator CMIteratorByIdType;

    CubitMeshSet_multiIndex &meshsets_index = m_mng->getMeshsetsMultindex();
    CMeshsetByType &m_by_type =
        meshsets_index.get<CubitMeshSets_mask_meshset_mi_tag>();
    CMeshsetByIdType &m_by_id_and_type =
        meshsets_index.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>();

    for (CMIteratorByType mit = m_by_type.lower_bound(SIDESET);
         mit != m_by_type.upper_bound(SIDESET); mit++) {

      std::cout << "Sideset on the mesh id = " << mit->getMeshsetId()
                << std::endl;
    }

    // iterate sideset and split
    for (int mm = 0; mm != nb_sidesets; mm++) {

      // find side set

      CMIteratorByIdType mit;
      mit = m_by_id_and_type.find(boost::make_tuple(sidesets[mm], SIDESET));
      if (mit == m_by_id_and_type.end()) {
        SETERRQ1(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
                 "No sideset in database id = %d", sidesets[mm]);
      }
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Split sideset %d\n",
                         mit->getMeshsetId());
      CHKERRQ(ierr);

      EntityHandle cubit_meshset = mit->getMeshset();
      {
        // get tet entities form back bit_level
        EntityHandle ref_level_meshset = 0;
        rval = moab.create_meshset(MESHSET_SET, ref_level_meshset);
        CHKERRQ_MOAB(rval);
        ierr = bit_mng->getEntitiesByTypeAndRefLevel(
            bit_levels.back(), BitRefLevel().set(), MBTET, ref_level_meshset);
        CHKERRQ(ierr);
        ierr = bit_mng->getEntitiesByTypeAndRefLevel(
            bit_levels.back(), BitRefLevel().set(), MBPRISM, ref_level_meshset);
        CHKERRQ(ierr);
        Range ref_level_tets;
        rval = moab.get_entities_by_handle(ref_level_meshset, ref_level_tets,
                                           true);
        CHKERRQ_MOAB(rval);

        // get faces and test to split
        ierr =
            interface_ptr->getSides(cubit_meshset, bit_levels.back(), true, 0);
        CHKERRQ(ierr);
        // set new bit level
        bit_levels.push_back(BitRefLevel().set(mm + 1));
        // split faces and
        ierr = interface_ptr->splitSides(ref_level_meshset, bit_levels.back(),
                                         cubit_meshset, false, true, 0);
        CHKERRQ(ierr);

        // clean meshsets
        rval = moab.delete_entities(&ref_level_meshset, 1);
        CHKERRQ_MOAB(rval);
      }
      // Update cubit meshsets
      for (_IT_CUBITMESHSETS_FOR_LOOP_(
        (core.getInterface<MeshsetsManager &, 0>()),ciit) ) {

        EntityHandle cubit_meshset = ciit->meshset;
        ierr = bit_mng->updateMeshsetByEntitiesChildren(
                       cubit_meshset, bit_levels.back(), cubit_meshset,
                       MBVERTEX, true);
        CHKERRQ(ierr);
        ierr = bit_mng->updateMeshsetByEntitiesChildren(
                       cubit_meshset, bit_levels.back(), cubit_meshset, MBEDGE,
                       true);
        CHKERRQ(ierr);
        ierr = bit_mng->updateMeshsetByEntitiesChildren(
                       cubit_meshset, bit_levels.back(), cubit_meshset, MBTRI,
                       true);
        CHKERRQ(ierr);
        ierr = bit_mng->updateMeshsetByEntitiesChildren(
                       cubit_meshset, bit_levels.back(), cubit_meshset, MBTET,
                       true);
        CHKERRQ(ierr);
        
      }
    }

    if (squash_bit_levels == PETSC_TRUE) {
      for (int ll = 0; ll != bit_levels.size() - 1; ll++) {
        ierr = m_field.delete_ents_by_bit_ref(bit_levels[ll], bit_levels[ll],
                                              true);
        CHKERRQ(ierr);
      }
      ierr = bit_mng->shiftRightBitRef(bit_levels.size() - 1);
      CHKERRQ(ierr);
    }

    if (output_vtk) {
      EntityHandle meshset;
      rval = moab.create_meshset(MESHSET_SET, meshset);
      CHKERRQ_MOAB(rval);
      BitRefLevel bit;
      if (squash_bit_levels)
        bit = bit_levels[0];
      else
        bit = bit_levels.back();
      ierr = bit_mng->getEntitiesByTypeAndRefLevel(
          bit, BitRefLevel().set(), MBTET, meshset);
      CHKERRQ(ierr);
      rval = moab.write_file("out.vtk", "VTK", "", &meshset, 1);
      CHKERRQ_MOAB(rval);
      rval = moab.delete_entities(&meshset, 1);
      CHKERRQ_MOAB(rval);
    }

    rval = moab.write_file("out.h5m");
    CHKERRQ_MOAB(rval);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF, e.errorCode, e.errorMessage);
  }

  ierr = PetscFinalize();
  CHKERRQ(ierr);

  return 0;
}
