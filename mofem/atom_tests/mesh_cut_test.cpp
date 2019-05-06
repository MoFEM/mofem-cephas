/** \file mesh_cut_test.cpp
 * \brief test for mesh cut interface
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

static char help[] = "testing mesh cut test\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);

    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "*** ERROR -my_file (MESH FILE NEEDED)");
    
    int side_set = 200;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-side_set", &side_set,
                              PETSC_NULL);

    double shift[] = {0, 0, 0};
    int nmax = 3;
    CHKERR PetscOptionsGetRealArray("", "-shift", shift, &nmax, &flg);
    if (flg && nmax != 3) 
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "three values expected");

    int nb_ref_before = 0;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-nb_ref_before", &nb_ref_before,
                              PETSC_NULL);
    int nb_ref_after = 0;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-nb_ref_after", &nb_ref_after,
                              PETSC_NULL);

    int fixed_edges_blockset = 100;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-fixed_edges_blockset",
                              &fixed_edges_blockset, PETSC_NULL);

    int corner_nodes_blockset = 1;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-corner_nodes_blockset",
                              &corner_nodes_blockset, PETSC_NULL);

    double tol[] = {0, 0, 0, 0};
    int nmax_tol = 4;
    CHKERR PetscOptionsGetRealArray("", "-tol", tol, &nmax_tol, &flg);
    if (flg && nmax_tol != 4) 
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "four values expected");

    PetscBool test = PETSC_FALSE;
    CHKERR
    PetscOptionsGetBool(PETSC_NULL, "", "-test", &test, PETSC_NULL);

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

    const char *option;
    option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    MoFEM::Core core(moab);
    MoFEM::CoreInterface &m_field =
        *(core.getInterface<MoFEM::CoreInterface>());

    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(
             (*core.getInterface<MeshsetsManager>()), SIDESET, it)) {
      cout << *it << endl;
    }

    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);
    BitRefLevel bit_last;
    bit_last.set(BITREFLEVEL_SIZE - 1);
    CHKERR m_field.getInterface<BitRefManager>()->addBitRefLevelByDim(0, 3,
                                                                      bit_last);
    CHKERR core.getInterface<BitRefManager>()->writeBitLevelByType(
        bit_last, BitRefLevel().set(), MBTET, "out_tets_bit_last.vtk", "VTK",
        "");

    int no_of_ents_not_in_database = -1;
    Range ents_not_in_database;
    if (test) {
      core.getInterface<BitRefManager>()->getAllEntitiesNotInDatabase(
          ents_not_in_database);
      no_of_ents_not_in_database = ents_not_in_database.size();
    }

    // get cut mesh interface
    CutMeshInterface *cut_mesh;
    CHKERR m_field.getInterface(cut_mesh);
    // get meshset manager interface
    MeshsetsManager *meshset_manager;
    CHKERR m_field.getInterface(meshset_manager);
    // get surface entities form side set
    Range surface;
    if (meshset_manager->checkMeshset(side_set, SIDESET)) 
      CHKERR meshset_manager->getEntitiesByDimension(side_set, SIDESET, 2,
                                                     surface, true);

    if (surface.empty()) 
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "No surface to cut");

    // Set surface entities. If surface entities are from existing side set,
    // copy those entities and do other geometrical transformations, like shift
    // scale or streach, rotate.
    if (meshset_manager->checkMeshset(side_set, SIDESET))
      CHKERR cut_mesh->copySurface(surface, NULL, shift, NULL, NULL,
                                   "surface.vtk");
    else
      CHKERR cut_mesh->setSurface(surface);

    // Get geometric corner nodes and corner edges
    Range fixed_edges, corner_nodes;
    if (meshset_manager->checkMeshset(fixed_edges_blockset, SIDESET)) {
      CHKERR meshset_manager->getEntitiesByDimension(
          fixed_edges_blockset, SIDESET, 1, fixed_edges, true);
    }
    if (meshset_manager->checkMeshset(fixed_edges_blockset, BLOCKSET)) {
      CHKERR meshset_manager->getEntitiesByDimension(
          fixed_edges_blockset, BLOCKSET, 1, fixed_edges, true);
    }
    if (meshset_manager->checkMeshset(corner_nodes_blockset, BLOCKSET)) {
      CHKERR meshset_manager->getEntitiesByDimension(
          corner_nodes_blockset, BLOCKSET, 0, corner_nodes, true);
    }

    // CHKERR cut_mesh->snapSurfaceToEdges(fixed_edges, 0.5, 0);

    Range tets;
    CHKERR moab.get_entities_by_dimension(0, 3, tets, false);
    CHKERR cut_mesh->setVolume(tets);

    // Build tree, it is used to ask geometrical queries, i.e. to find edges to
    // cut or trim.
    CHKERR cut_mesh->buildTree();

    // Create tag storing nodal positions
    double def_position[] = {0, 0, 0};
    Tag th;
    CHKERR moab.tag_get_handle("POSITION", 3, MB_TYPE_DOUBLE, th,
                               MB_TAG_CREAT | MB_TAG_SPARSE, def_position);
    // Set tag values with coordinates of nodes
    CHKERR cut_mesh->setTagData(th);

    // Get BitRefManager interface,,
    BitRefManager *bit_ref_manager;
    CHKERR m_field.getInterface(bit_ref_manager);

    // Cut mesh, trim surface and merge bad edges
    int first_bit = 1;
    CHKERR cut_mesh->refCutTrimAndMerge(
        first_bit, 1, nb_ref_before, nb_ref_after, th, tol[0], tol[1], tol[2],
        tol[3], fixed_edges, corner_nodes, true, true);

    if (test) {
      struct TestBitLevel {
        BitRefManager *mngPtr;
        TestBitLevel(BitRefManager *mng_ptr) : mngPtr(mng_ptr) {}
        MoFEMErrorCode operator()(const BitRefLevel &bit,
                                  const int expected_size) {
          MoFEMFunctionBeginHot;
          Range ents;
          CHKERR mngPtr->getEntitiesByRefLevel(bit, bit, ents);
          cout << "bit_level nb ents " << ents.size() << endl;
          if (expected_size != -1 &&
              expected_size != static_cast<int>(ents.size())) {
            SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                     "Wrong bit ref size %d!=%d", expected_size, ents.size());
          }
          MoFEMFunctionReturnHot(0);
        }
      };
      for (int ll = 1; ll != first_bit; ++ll)
        CHKERR TestBitLevel(core.getInterface<BitRefManager>())(
            BitRefLevel().set(ll), -1);
    }

    // Improve mesh with tetgen
#undef WITH_TETGEN
#ifdef WITH_TETGEN
    int bit_tetgen = first_bit + 1;
    // Switches controling TetGen
    vector<string> switches;
    switches.push_back("YrqOJMS0VV");
    CHKERR cut_mesh->rebuildMeshWithTetGen(
        switches, BitRefLevel().set(first_bit), BitRefLevel().set(bit_tetgen),
        cut_mesh->getMergedSurfaces(), fixed_edges, corner_nodes, th, true);
    CHKERR core.getInterface<BitRefManager>()->writeBitLevelByType(
        BitRefLevel().set(bit_tetgen), BitRefLevel().set(), MBTET,
        "out_tets_tetgen.vtk", "VTK", "");
#else
    int bit_tetgen = first_bit;
    const_cast<Range &>(cut_mesh->getTetgenSurfaces()) =
        cut_mesh->getMergedSurfaces();
#endif // WITH_TETGEN

    // Split faces
    CHKERR cut_mesh->splitSides(BitRefLevel().set(bit_tetgen),
                                BitRefLevel().set(bit_tetgen + 1),
                                cut_mesh->getTetgenSurfaces(), th);
    CHKERR core.getInterface<MeshsetsManager>()
        ->updateAllMeshsetsByEntitiesChildren(
            BitRefLevel().set(bit_tetgen + 1));

    CHKERR core.getInterface<BitRefManager>()->writeBitLevelByType(
        BitRefLevel().set(bit_tetgen + 1), BitRefLevel().set(), MBTET,
        "out_split.vtk", "VTK", "");

    // Finally shift bits
    BitRefLevel shift_mask;
    for (int ll = 0; ll != bit_tetgen + 2; ++ll)
      shift_mask.set(ll);
    CHKERR core.getInterface<BitRefManager>()->shiftRightBitRef(
        bit_tetgen, shift_mask, VERBOSE);

    // Set coordinates for tag data
    CHKERR cut_mesh->setCoords(th);
    CHKERR core.getInterface<BitRefManager>()->writeBitLevelByType(
        BitRefLevel().set(0), BitRefLevel().set(), MBTET,
        "out_tets_shift_level0.vtk", "VTK", "");
    CHKERR core.getInterface<BitRefManager>()->writeBitLevelByType(
        BitRefLevel().set(1), BitRefLevel().set(), MBTET,
        "out_tets_shift_level1.vtk", "VTK", "");

    Range surface_verts;
    CHKERR moab.get_connectivity(cut_mesh->getSurface(), surface_verts);
    Range adj_surface_edges;
    CHKERR moab.get_adjacencies(cut_mesh->getSurface(), 1, false,
                                adj_surface_edges, moab::Interface::UNION);
    CHKERR moab.delete_entities(cut_mesh->getSurface());
    CHKERR moab.delete_entities(adj_surface_edges);
    CHKERR moab.delete_entities(surface_verts);
    CHKERR m_field.delete_ents_by_bit_ref(
        BitRefLevel().set(0) | BitRefLevel().set(1),
        BitRefLevel().set(0) | BitRefLevel().set(1), true, VERBOSE);

    {
      EntityHandle meshset;
      CHKERR moab.create_meshset(MESHSET_SET, meshset);
      Range tets;
      CHKERR moab.get_entities_by_dimension(0, 3, tets, true);
      CHKERR moab.add_entities(meshset, tets);
      CHKERR moab.write_file("out.vtk", "VTK", "", &meshset, 1);
      CHKERR moab.delete_entities(&meshset, 1);
    }
    CHKERR core.getInterface<BitRefManager>()->writeBitLevelByType(
        bit_last, BitRefLevel().set(), MBTET, "out_tets_bit_last.vtk", "VTK",
        "");
    CHKERR core.getInterface<BitRefManager>()->writeEntitiesNotInDatabase(
        "left_entities.vtk", "VTK", "");

    if (test) {
      Range ents;
      core.getInterface<BitRefManager>()->getAllEntitiesNotInDatabase(ents);
      if (no_of_ents_not_in_database != static_cast<int>(ents.size())) {
        cerr << subtract(ents, ents_not_in_database) << endl;
        EntityHandle meshset;
        CHKERR moab.create_meshset(MESHSET_SET, meshset);
        Range tets;
        CHKERR moab.get_entities_by_dimension(0, 3, tets, true);
        CHKERR moab.add_entities(meshset, subtract(ents, ents_not_in_database));
        CHKERR moab.write_file("not_cleanded.vtk", "VTK", "", &meshset, 1);
        CHKERR moab.delete_entities(&meshset, 1);
        SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                 "Inconsistent number of ents %d!=%d",
                 no_of_ents_not_in_database, ents.size());
      }
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
