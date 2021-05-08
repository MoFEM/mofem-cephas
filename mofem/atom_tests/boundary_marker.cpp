/**
 * \file boundary_marker.cpp
 * \example boundary_marker.cpp
 *
 * Mark skin entities, i.e. boundary, next check DOFs if are properly marked
 * when accessed form user data operator.
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

using FaceEle = MoFEM::FaceElementForcesAndSourcesCore;
using FaceEleOp = FaceEle::UserDataOperator;

struct OpFace : public FaceEleOp {
  const Range &skinEnts;
  const std::vector<unsigned char> &mArker;

  OpFace(const Range &skin_ents, const std::vector<unsigned char> &marker)
      : FaceEleOp("FIELD1", OPROW), skinEnts(skin_ents), mArker(marker) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;

    const int nb_dofs = data.getIndices().size();
    if (nb_dofs == 0)
      MoFEMFunctionReturnHot(0);

    // This function takes entities on DOFs, and check if those entities are
    // contained in the range. Note it is local element vector, which is used
    // to validate of global local vector.
    auto get_boundary_marker_directly_from_range = [&]() {
      std::vector<unsigned char> ents_marker_used_for_testing;
      ents_marker_used_for_testing.resize(data.getLocalIndices().size(),0);
      switch (type) {
      case MBVERTEX: {
        for (size_t side = 0; side != getNumberOfNodesOnElement(); ++side) {
          EntityHandle ent = getSideEntity(side, MBVERTEX);
          ents_marker_used_for_testing[3 * side + 0] =
              skinEnts.find(ent) != skinEnts.end() ? 1 : 0;
          ents_marker_used_for_testing[3 * side + 2] =
              ents_marker_used_for_testing[3 * side + 0];
        }
        break;
      }
      default: {
        EntityHandle ent = getSideEntity(side, type);
        bool is_on_boundary = skinEnts.find(ent) != skinEnts.end();
        for (size_t dof = 0; dof != nb_dofs; ++dof)
          if ((dof % 3) != 1)
            ents_marker_used_for_testing[dof] = is_on_boundary ? 1 : 0;
      }
      };
      return ents_marker_used_for_testing;
    };

    auto test_marker = get_boundary_marker_directly_from_range();
    for (size_t n = 0; n != nb_dofs; ++n) {
      const size_t local_index = data.getLocalIndices()[n];
      if (test_marker[n] != mArker[local_index])
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "Wrong boundary marker");
    }

    MoFEMFunctionReturn(0);
  }
};

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

    int order = 4;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);

    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    auto *simple_interface = m_field.getInterface<Simple>();
    CHKERR simple_interface->getOptions();
    CHKERR simple_interface->loadFile("");
    CHKERR simple_interface->addDomainField("FIELD1", H1,
                                            AINSWORTH_LEGENDRE_BASE, 3);
    CHKERR simple_interface->setFieldOrder("FIELD1", order);
    CHKERR simple_interface->setUp();

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

    auto mark_boundary_dofs = [&](Range &skin_edges) {
      auto problem_manager = m_field.getInterface<ProblemsManager>();
      std::vector<unsigned char> marker;
      problem_manager->markDofs(simple_interface->getProblemName(), ROW,
                                ProblemsManager::OR, skin_edges, marker);
      // Unset coefficient 1, e.g. u_y
      problem_manager->modifyMarkDofs(simple_interface->getProblemName(), ROW,
                                      "FIELD1", 1, 1, ProblemsManager::AND, 0,
                                      marker);
      return marker;
    };

    auto skin_ents = get_ents_on_mesh_skin();

    // Get global local vector of marked DOFs. Is global, since is set for all
    // DOFs on processor. Is local since only DOFs on processor are in the
    // vector. To access DOFs use local indices.
    auto marker = mark_boundary_dofs(skin_ents);

    boost::shared_ptr<FaceEle> fe(new FaceEle(m_field));
    fe->getOpPtrVector().push_back(new OpFace(skin_ents, marker));

    auto dm = simple_interface->getDM();
    CHKERR DMoFEMLoopFiniteElements(dm, simple_interface->getDomainFEName(),
                                    fe);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
