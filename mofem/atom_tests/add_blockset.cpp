/**
 * \file add_blockset.cpp
 * \example add_blockset.cpp
 *
 * Create blockset and add entities. Next check if entities are in the blockset.
 *
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

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

    // Create MoFEM instance
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;
    const char *option = "";
    CHKERR moab.load_file("rectangle.h5m", 0, option);

    auto meshsets_mng = m_field.getInterface<MeshsetsManager>();

    auto get_ents_on_mesh_skin = [&]() {
      Range faces;
      CHKERR m_field.get_moab().get_entities_by_type(0, MBTRI, faces);
      Skinner skin(&m_field.get_moab());
      Range skin_edges;
      CHKERR skin.find_skin(0, faces, false, skin_edges);
      Range skin_verts;
      CHKERR moab.get_connectivity(skin_edges, skin_verts, true);
      skin_edges.merge(skin_verts);
      return skin_edges;
    };

    auto add_blockset = [&](const Range skin_ents) {
      MoFEMFunctionBegin;
      CHKERR meshsets_mng->addMeshset(BLOCKSET, 1);
      CHKERR meshsets_mng->addEntitiesToMeshset(BLOCKSET, 1, skin_ents);
      MoFEMFunctionReturn(0);
    };

    auto print_blocksets = [&]() {
      MoFEMFunctionBegin;
      for (auto &it : meshsets_mng->getMeshsetsMultindex())
        cout << it << endl;
      MoFEMFunctionReturn(0);
    };

    auto check_meshset = [&](const Range skin_ents) {
      MoFEMFunctionBegin;
      std::vector<EntityHandle> ents(skin_ents.size());
      std::copy(skin_ents.begin(), skin_ents.end(), ents.begin());
      const bool test = meshsets_mng->checkIfMeshsetContainsEntities(
          1, BLOCKSET, &*ents.begin(), ents.size());
      if (!test)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "All entities should be in blockset");
      MoFEMFunctionReturn(0);
    };

    auto skin_ents = get_ents_on_mesh_skin();
    CHKERR add_blockset(skin_ents);
    CHKERR print_blocksets();
    CHKERR check_meshset(skin_ents);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
