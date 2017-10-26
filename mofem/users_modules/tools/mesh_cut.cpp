/** \file mesh_cut.cpp
 * \brief test for mesh cut interface
 * \example mesh_cut.cpp
 *
 * \ingroup mesh_cut
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

static char help[] = "mesh cutting\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc, &argv, (char *)0, help);

  try {

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
    int surface_side_set = 200;
    PetscBool flg_vol_block_set;
    int vol_block_set = 1;
    int edges_block_set = 1;
    int vertex_block_set = 2;
    double shift[] = {0, 0, 0};
    int nmax = 3;
    int fraction_level = 2;
    PetscBool squash_bits = PETSC_TRUE;
    PetscBool set_coords = PETSC_TRUE;

    ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "", "Mesh cut options", "none");
    CHKERRQ(ierr);

    ierr = PetscOptionsString("-my_file", "mesh file name", "", "mesh.h5m",
                              mesh_file_name, 255, &flg);
    CHKERRQ(ierr);
    ierr = PetscOptionsInt("-surface_side_set", "surface side set", "",
                           surface_side_set, &surface_side_set, PETSC_NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsInt("-vol_block_set", "volume side set", "",
                           vol_block_set, &vol_block_set, &flg_vol_block_set);
    CHKERRQ(ierr);
    ierr = PetscOptionsInt("-edges_block_set", "edges side set", "",
                           edges_block_set, &edges_block_set, PETSC_NULL);
    CHKERRQ(ierr);

    ierr = PetscOptionsInt("-vertex_block_set", "vertex side set", "",
                           vertex_block_set, &vertex_block_set, PETSC_NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsRealArray("-shift", "shift surface by vector", "", shift,
                                 &nmax, &flg);
    CHKERRQ(ierr);
    ierr = PetscOptionsInt("-fraction_level", "fraction of merges merged", "",
                           fraction_level, &fraction_level, PETSC_NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsBool("-squash_bits", "true to squash bits at the end",
                            "", squash_bits, &squash_bits, PETSC_NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsBool("-set_coords", "true to set coords at the end", "",
                            set_coords, &set_coords, PETSC_NULL);
    CHKERRQ(ierr);

    ierr = PetscOptionsEnd();
    CHKERRQ(ierr);

    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "*** ERROR -my_file (MESH FILE NEEDED)");
    }
    if (flg && nmax != 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "three values expected");
    }
 

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

    const char *option;
    option = ""; //"PARALLEL=BCAST";//;DEBUG_IO";
    rval = moab.load_file(mesh_file_name, 0, option);
    CHKERRQ_MOAB(rval);

    MoFEM::Core core(moab);
    MoFEM::CoreInterface &m_field =
        *(core.getInterface<MoFEM::CoreInterface>());

      // get cut mesh interface
    CutMeshInterface *cut_mesh;
    ierr = m_field.getInterface(cut_mesh);
    CHKERRQ(ierr);
    // get meshset manager interface
    MeshsetsManager *meshset_manager;
    ierr = m_field.getInterface(meshset_manager);
    CHKERRQ(ierr);
    // get bit ref manager interface
    BitRefManager *bit_ref_manager;
    ierr = m_field.getInterface(bit_ref_manager); CHKERRQ(ierr);

    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(
             (core.getInterface<MeshsetsManager &, 0>()), SIDESET, it)) {
      cout << *it << endl;
    }

    BitRefLevel bit_level0;
    bit_level0.set(0);
    ierr = bit_ref_manager->setBitRefLevelByDim(0, 3, bit_level0);
    CHKERRQ(ierr);

    // get surface entities form side set
    Range surface;
    if (meshset_manager->checkMeshset(surface_side_set, SIDESET)) {
      ierr = meshset_manager->getEntitiesByDimension(surface_side_set, SIDESET, 2,
                                                     surface, true);
      CHKERRQ(ierr);
    }
    if (surface.empty()) {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "No surface to cut");
    }

    // Set surface entities. If surface entities are from existing side set,
    // copy those entities and do other geometrical transformations, like shift
    // scale or streach, rotate.
    if (meshset_manager->checkMeshset(surface_side_set, SIDESET)) {
      ierr = cut_mesh->copySurface(surface, NULL, shift);
      CHKERRQ(ierr);
    } else {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Side set not found %d", surface_side_set);
    }

    Range tets;
    if (flg_vol_block_set) {
      if (meshset_manager->checkMeshset(vol_block_set, BLOCKSET)) {
        ierr = meshset_manager->getEntitiesByDimension(
            vol_block_set, BLOCKSET, 3, tets, true);
        CHKERRQ(ierr);
      } else {
        SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Block set %d not found", vol_block_set);
      }
    } else {
      ierr = moab.get_entities_by_dimension(0, 3, tets, false);
      CHKERRQ(ierr);
    }
    ierr = cut_mesh->setVolume(tets);
    CHKERRQ(ierr);

    // Build tree, it is used to ask geometrical queries, i.e. to find edges to
    // cut or trim.
    ierr = cut_mesh->buildTree();
    CHKERRQ(ierr);

    BitRefLevel bit_level1; // Cut level
    bit_level1.set(1);
    BitRefLevel bit_level2; // Trim level
    bit_level2.set(2);
    BitRefLevel bit_level3; // Merge level
    bit_level3.set(3);
    BitRefLevel bit_level4; // TeteGen level
    bit_level4.set(4);

    // Create tag storing nodal positions
    double def_position[] = {0, 0, 0};
    Tag th;
    rval = moab.tag_get_handle("POSITION", 3, MB_TYPE_DOUBLE, th,
                               MB_TAG_CREAT | MB_TAG_SPARSE, def_position);
    CHKERRQ_MOAB(rval);
    // Set tag values with coordinates of nodes
    ierr = cut_mesh->setTagData(th);
    CHKERRQ(ierr);

    // Get geometric corner nodes and corner edges
    Range fixed_edges, corner_nodes;
    if (meshset_manager->checkMeshset(edges_block_set, BLOCKSET)) {
      ierr = meshset_manager->getEntitiesByDimension(edges_block_set, BLOCKSET,
                                                     1, fixed_edges, true);
      CHKERRQ(ierr);
    }
    if (meshset_manager->checkMeshset(vertex_block_set, BLOCKSET)) {
      ierr = meshset_manager->getEntitiesByDimension(
          vertex_block_set, BLOCKSET, 0, corner_nodes, true);
      CHKERRQ(ierr);
    }

    // Cut mesh, trim surface and merge bad edges
    ierr = cut_mesh->cutTrimAndMerge(fraction_level, bit_level1, bit_level2,
                                     bit_level3, th, 1e-4, 1e-2, 1e-4, 1e-2,
                                     fixed_edges, corner_nodes, true, false);
    CHKERRQ(ierr);

    // Improve mesh with tetgen
#ifdef WITH_TETGEN
    // Switches controling TetGen
    vector<string> switches;
    switches.push_back("rp180YsqORJS0VV");
    ierr = cut_mesh->rebuildMeshWithTetGen(switches, bit_level3, bit_level4,
                                           cut_mesh->getMergedSurfaces(),
                                           fixed_edges, corner_nodes, th);
    CHKERRQ(ierr);
#else
    bit_level4 = bit_level3;
    const_cast<Range &>(cut_mesh->getTetgenSurfaces()) =
        cut_mesh->getMergedSurfaces();
#endif // WITH_TETGEN

    // Add tets from last level to block
    if(flg_vol_block_set) {
      EntityHandle meshset;
      ierr = meshset_manager->getMeshset(vol_block_set, BLOCKSET, meshset);
      CHKERRQ(ierr);
      ierr = bit_ref_manager->getEntitiesByTypeAndRefLevel(
          bit_level4, BitRefLevel().set(), MBTET, meshset);
      CHKERRQ(ierr);
    }

    // Shift bits
    if (squash_bits) {
      BitRefLevel shift_mask;
      shift_mask.set(0);
      shift_mask.set(1);
      shift_mask.set(2);
      shift_mask.set(3);
      shift_mask.set(4);
      ierr = core.getInterface<BitRefManager>()->shiftRightBitRef(4, shift_mask,
                                                                  VERBOSE);
      CHKERRQ(ierr);
    }

    if (set_coords) {
      // Set coordinates for tag data
      ierr = cut_mesh->setCoords(th);
      CHKERRQ(ierr);
    }

    Range surface_verts;
    rval = moab.get_connectivity(cut_mesh->getSurface(), surface_verts);
    CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(cut_mesh->getSurface());
    CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(surface_verts);
    CHKERRQ_MOAB(rval);

    rval = moab.write_file("out.h5m"); CHKERRQ_MOAB(rval);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF, e.errorCode, e.errorMessage);
  }

  PetscFinalize();

  return 0;
}