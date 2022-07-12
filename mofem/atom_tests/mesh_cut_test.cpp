/** \file mesh_cut_test.cpp
 * \brief test for mesh cut interface
 *
 * \ingroup mesh_cut
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "testing mesh cut test\n\n";

struct TestBitLevel {
  BitRefManager *mngPtr;
  TestBitLevel(BitRefManager *mng_ptr) : mngPtr(mng_ptr) {}
  MoFEMErrorCode operator()(const BitRefLevel &bit, const int expected_size) {
    MoFEMFunctionBeginHot;
    Range ents;
    CHKERR mngPtr->getEntitiesByRefLevel(bit, BitRefLevel().set(), ents);
    MOFEM_LOG("WORLD", Sev::inform)
        << "bit_level nb ents " << bit << " " << ents.size();
    if (expected_size != -1 && expected_size != static_cast<int>(ents.size())) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
               "Wrong bit ref size %d!=%d", expected_size, ents.size());
    }
    MoFEMFunctionReturnHot(0);
  }
};

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    MOFEM_LOG_CHANNEL("WORLD");

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

    int restricted_side_set = 205;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-restricted_side_set",
                              &restricted_side_set, PETSC_NULL);

    double shift[] = {0, 0, 0};
    int nmax = 3;
    CHKERR PetscOptionsGetRealArray("", "-shift", shift, &nmax, &flg);
    if (flg && nmax != 3)
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "three values expected");

    int fixed_edges_blockset = 100;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-fixed_edges_blockset",
                              &fixed_edges_blockset, PETSC_NULL);

    int corner_nodes_blockset = 1;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-corner_nodes_blockset",
                              &corner_nodes_blockset, PETSC_NULL);

    double tol[] = {0, 0, 0};
    int nmax_tol = 3;
    CHKERR PetscOptionsGetRealArray("", "-tol", tol, &nmax_tol, &flg);
    if (flg && nmax_tol != 3)
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "four values expected");

    PetscBool test = PETSC_FALSE;
    CHKERR
    PetscOptionsGetBool(PETSC_NULL, "", "-test", &test, PETSC_NULL);

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    auto moab_comm_wrap =
        boost::make_shared<WrapMPIComm>(PETSC_COMM_WORLD, false);
    if (pcomm == NULL)
      pcomm =
          new ParallelComm(&moab, moab_comm_wrap->get_comm());

    const char *option;
    option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    MoFEM::Core core(moab);
    MoFEM::CoreInterface &m_field =
        *(core.getInterface<MoFEM::CoreInterface>());

    MOFEM_LOG_CHANNEL("WORLD");
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(
             (*core.getInterface<MeshsetsManager>()), SIDESET, it)) {
      MOFEM_LOG("WORLD", Sev::inform) << *it;
    }

    BitRefLevel bit_level0;
    bit_level0.set(0);
    BitRefLevel bit_last;
    bit_last.set(BITREFLEVEL_SIZE - 1);

    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, BitRefLevel().set(0));
    CHKERR m_field.getInterface<BitRefManager>()->addBitRefLevelByDim(0, 3,
                                                                      bit_last);

    CHKERR core.getInterface<BitRefManager>()->writeBitLevelByType(
        bit_level0, BitRefLevel().set(), MBTET, "out_tets_init_level0.vtk",
        "VTK", "");
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

    // Get BitRefManager interface,,
    BitRefManager *bit_ref_manager;
    CHKERR m_field.getInterface(bit_ref_manager);
    // get meshset manager interface
    MeshsetsManager *meshset_manager;
    CHKERR m_field.getInterface(meshset_manager);
    // get cut mesh interface
    CutMeshInterface *cut_mesh;
    CHKERR m_field.getInterface(cut_mesh);

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

    // get surface entities form side set
    Range restriced_surface;
    if (meshset_manager->checkMeshset(restricted_side_set, SIDESET))
      CHKERR meshset_manager->getEntitiesByDimension(
          restricted_side_set, SIDESET, 2, restriced_surface, true);

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

    cut_mesh->setConstrainSurface(restriced_surface);

    Range tets;
    CHKERR moab.get_entities_by_dimension(0, 3, tets, false);

    CHKERR cut_mesh->setVolume(tets);
    CHKERR cut_mesh->buildTree();
    CHKERR cut_mesh->makeFront(true);
    const int nb_ref_cut = 1;
    const int nb_ref_trim = 1;
    CHKERR cut_mesh->refineMesh(0, nb_ref_cut, nb_ref_trim, &fixed_edges,
                                VERBOSE, false);
    auto shift_after_ref = [&]() {
      MoFEMFunctionBegin;
      BitRefLevel mask;
      mask.set(0);
      for (int ll = 1; ll != nb_ref_cut + nb_ref_trim + 1; ++ll)
        mask.set(ll);
      CHKERR core.getInterface<BitRefManager>()->shiftRightBitRef(
          nb_ref_cut + nb_ref_trim, mask, VERBOSE);
      MoFEMFunctionReturn(0);
    };
    CHKERR shift_after_ref();

    CHKERR m_field.getInterface<BitRefManager>()
        ->writeEntitiesAllBitLevelsByType(BitRefLevel().set(), MBTET,
                                          "all_bits.vtk", "VTK", "");

    // Create tag storing nodal positions
    double def_position[] = {0, 0, 0};
    Tag th;
    CHKERR moab.tag_get_handle("POSITION", 3, MB_TYPE_DOUBLE, th,
                               MB_TAG_CREAT | MB_TAG_SPARSE, def_position);
    // Set tag values with coordinates of nodes
    CHKERR cut_mesh->setTagData(th);

    // Cut mesh, trim surface and merge bad edges
    int first_bit = 1;
    CHKERR cut_mesh->cutTrimAndMerge(first_bit, 5, th, tol[0], tol[1], tol[2],
                                     fixed_edges, corner_nodes, true, false);

    // Split faces
    CHKERR cut_mesh->splitSides(BitRefLevel().set(first_bit),
                                BitRefLevel().set(first_bit + 1),
                                cut_mesh->getMergedSurfaces(), th);
    CHKERR core.getInterface<MeshsetsManager>()
        ->updateAllMeshsetsByEntitiesChildren(BitRefLevel().set(first_bit + 1));

    CHKERR core.getInterface<BitRefManager>()->writeBitLevelByType(
        BitRefLevel().set(first_bit + 1), BitRefLevel().set(), MBTET,
        "out_split_tets.vtk", "VTK", "");

    CHKERR core.getInterface<BitRefManager>()->writeBitLevelByType(
        BitRefLevel().set(first_bit + 1), BitRefLevel().set(), MBPRISM,
        "out_split_prism.vtk", "VTK", "");

    if (test) {
      for (int ll = 0; ll != first_bit + 2; ++ll)
        CHKERR TestBitLevel(core.getInterface<BitRefManager>())(
            BitRefLevel().set(ll), -1);
    }

    // Finally shift bits
    BitRefLevel shift_mask;
    for (int ll = 0; ll != first_bit; ++ll)
      shift_mask.set(ll);
    CHKERR core.getInterface<BitRefManager>()->shiftRightBitRef(
        first_bit, shift_mask, VERBOSE);

    // Set coordinates for tag data
    CHKERR cut_mesh->setCoords(th);
    CHKERR core.getInterface<BitRefManager>()->writeBitLevelByType(
        BitRefLevel().set(first_bit), BitRefLevel().set(), MBTET,
        "out_tets_shift_level0.vtk", "VTK", "");
    CHKERR core.getInterface<BitRefManager>()->writeBitLevelByType(
        BitRefLevel().set(first_bit + 1), BitRefLevel().set(), MBTET,
        "out_tets_shift_level1.vtk", "VTK", "");
    CHKERR core.getInterface<BitRefManager>()->writeBitLevelByType(
        BitRefLevel().set(first_bit + 1), BitRefLevel().set(), MBPRISM,
        "out_tets_shift_level1_prism.vtk", "VTK", "");

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

    MOFEM_LOG_CHANNEL("WORLD");
    if (test) {
      Range ents;
      core.getInterface<BitRefManager>()->getAllEntitiesNotInDatabase(ents);
      if (no_of_ents_not_in_database != static_cast<int>(ents.size())) {
        MOFEM_LOG("WORLD", Sev::inform) << subtract(ents, ents_not_in_database);
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
