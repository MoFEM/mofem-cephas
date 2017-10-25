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
    ierr = PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
    CHKERRQ(ierr);
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "*** ERROR -my_file (MESH FILE NEEDED)");
    }
    int surface_side_set = 200;
    ierr =
        PetscOptionsGetInt(PETSC_NULL, "", "-surface_side_set", &surface_side_set, PETSC_NULL);
    CHKERRQ(ierr);
    int edges_block_set = 1;
    ierr = PetscOptionsGetInt(PETSC_NULL, "", "-edges_block_set", &edges_block_set,
                              PETSC_NULL);
    CHKERRQ(ierr);
    int vertex_block_set = 2;
    ierr = PetscOptionsGetInt(PETSC_NULL, "", "-vertex_block_set", &vertex_block_set,
                              PETSC_NULL);

    double shift[] = {0, 0, 0};
    int nmax = 3;
    ierr = PetscOptionsGetRealArray("", "-shift", shift, &nmax, &flg);
    if (flg && nmax != 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "three values expected");
    }
    int fraction_level = 2;
    ierr =
        PetscOptionsGetInt(PETSC_NULL, "", "-fraction_level", &fraction_level, PETSC_NULL);
    CHKERRQ(ierr);
    PetscBool squash_bits = PETSC_TRUE;
    ierr =
        PetscOptionsGetBool(PETSC_NULL, "", "-squash_bits", &squash_bits, PETSC_NULL);
    CHKERRQ(ierr);

    PetscBool set_coords = PETSC_TRUE;
    ierr =
        PetscOptionsGetBool(PETSC_NULL, "", "-set_coords", &set_coords, PETSC_NULL);
    CHKERRQ(ierr);

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

    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(
             (*core.getInterface<MeshsetsManager>()), SIDESET, it)) {
      cout << *it << endl;
    }

    BitRefLevel bit_level0;
    bit_level0.set(0);
    ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);
    CHKERRQ(ierr);

    // get cut mesh interface
    CutMeshInterface *cut_mesh;
    ierr = m_field.getInterface(cut_mesh);
    CHKERRQ(ierr);

    // get meshset manager interface
    MeshsetsManager *meshset_manager;
    ierr = m_field.getInterface(meshset_manager);
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
    } 

    Range tets;
    ierr = moab.get_entities_by_dimension(0, 3, tets, false);
    CHKERRQ(ierr);
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

    // Get BitRefManager interface
    BitRefManager *bit_ref_manager;
    ierr = m_field.getInterface(bit_ref_manager);
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