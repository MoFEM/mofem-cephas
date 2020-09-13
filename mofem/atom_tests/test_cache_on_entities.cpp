/**
 * \file test_cache_entities.cpp
 * \example test_cache_entities.cpp
 *
 * Tetsing entities cache. Entities acache is data abstraction enabling user to
 * pass data between operators, finite elements, and problems, in convinient
 * way.
 *
 * It can be used to store indices, or history variables, and any other data,
 * which I can not imagine at that moment of time, and how to make use of it
 * depends on user creativity and imagination.
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

using VolEle = VolumeElementForcesAndSourcesCore;
using VolOp = VolEle::UserDataOperator;
using FaceEle = FaceElementForcesAndSourcesCore;
using FaceOp = FaceEle::UserDataOperator;

struct MyStorage : EntityStorage {
  MyStorage(VectorInt &data) : globalIndices(data) {}
  VectorInt globalIndices;
};

/**
 * @brief Operator set cache stored data, in this example, global indices, but
 * it can be any structure
 *
 */
struct OpVolumeSet : public VolOp {
  OpVolumeSet(const std::string &field_name) : VolOp(field_name, OPROW) {}
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {

    MoFEMFunctionBegin;
    MOFEM_LOG_CHANNEL("SYNC");
    MOFEM_LOG_TAG("SYNC", "OpVolumeSet");

    // Clear data when start process element
    if (type == MBVERTEX)
      entsIndices.clear();

    // Get field entity pointer
    auto field_ents = data.getFieldEntities();
    if (auto e_ptr = field_ents[0]) {
      // Add indices to global storage
      entsIndices.push_back(boost::make_shared<MyStorage>(data.getIndices()));
      // Store pointer to data on entity
      e_ptr->getWeakStoragePtr() = entsIndices.back();

      MOFEM_LOG("SYNC", Sev::inform) << "Set " << e_ptr->getEntTypeName() << " "
                                     << side << " : " << entsIndices.size();

      // Check if all works
      if (auto ptr = e_ptr->getSharedStoragePtr<EntityStorage>()) {

        if (auto cast_ptr = boost::dynamic_pointer_cast<MyStorage>(ptr))
          MOFEM_LOG("SYNC", Sev::noisy) << "Cast works";
        else
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Cast not works");

      } else {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Pointer to stored data is expected to be set");
      }

      if (auto ptr = e_ptr->getSharedStoragePtr<MyStorage>()) {
        MOFEM_LOG("SYNC", Sev::verbose)
            << data.getIndices() << " : " << ptr->globalIndices;
      } else {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Pointer to stored data is expected to be set");
      }

    } else {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Pointer to entity should be set by finite element");
    }

    MoFEMFunctionReturn(0);
  }

  using SharedVecInt = boost::shared_ptr<MyStorage>;
  static std::vector<SharedVecInt>
      entsIndices; ///< This is global static storage
};

std::vector<OpVolumeSet::SharedVecInt> OpVolumeSet::entsIndices;

/**
 * @brief Test if cached data can be accessed, and check consistency of data
 *
 */
struct OpVolumeTest : public VolOp {
  OpVolumeTest(const std::string &field_name) : VolOp(field_name, OPROW) {}
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {

    MoFEMFunctionBegin;
    MOFEM_LOG_CHANNEL("SYNC");
    MOFEM_LOG_TAG("SYNC", "OpVolumeTest");

    // Get pointer to field entities
    auto field_ents = data.getFieldEntities();
    if (auto e_ptr = field_ents[0]) {

      MOFEM_LOG("SYNC", Sev::inform)
          << "Test " << e_ptr->getEntTypeName() << " " << side;

      // Check if data are cached on entity, and if code is correct, data should
      // accessible.
      if (auto ptr = e_ptr->getSharedStoragePtr<MyStorage>()) {

        MOFEM_LOG("SYNC", Sev::verbose)
            << data.getIndices() << " : " << ptr->globalIndices;

        // Check constancy of data. Stored data are indices, and expected stored
        // that should be indices, thus difference between should be zero.
        auto diff = data.getIndices() - ptr->globalIndices;
        auto dot = inner_prod(diff, diff);
        if (dot > 0)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Test failed");

      } else {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Pointer to stored data is expected to be set");
      }

    } else {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Pointer to entity should be set by finite element");
    }

    MoFEMFunctionReturn(0);
  }
};

struct OpVolumeAssemble : public VolOp {
  OpVolumeAssemble(const std::string &field_name, Vec v)
      : VolOp(field_name, OPROW), V(v) {}
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    auto nb_dofs = data.getIndices().size();
    if (nb_dofs) {
      std::vector<double> v(nb_dofs, 1);
      CHKERR VecSetValues<EssentialBcStorage>(V, data, &v[0], ADD_VALUES);
    }
    MoFEMFunctionReturn(0);
  }

private:
  Vec V;
};

struct VolRule {
  int operator()(int, int, int) const { return 1; }
};

int main(int argc, char *argv[]) {

  // initialize petsc
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    // Create MoAB database
    moab::Core moab_core;
    moab::Interface &moab = moab_core;

    // Create MoFEM database and link it to MoAB
    MoFEM::Core mofem_core(moab);
    MoFEM::Interface &m_field = mofem_core;

    // Register DM Manager
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Simple interface
    Simple *simple_interface;
    CHKERR m_field.getInterface(simple_interface);
    {
      // get options from command line
      CHKERR simple_interface->getOptions();
      // load mesh file
      CHKERR simple_interface->loadFile();
      // add fields
      CHKERR simple_interface->addDomainField("FIELD", H1,
                                              AINSWORTH_LEGENDRE_BASE, 3);
      // set fields order
      CHKERR simple_interface->setFieldOrder("FIELD", 4);
      // setup problem
      CHKERR simple_interface->setUp();

      // get dm
      auto dm = simple_interface->getDM();
      // create elements
      auto domain_fe = boost::make_shared<VolEle>(m_field);
      // set integration rule
      domain_fe->getRuleHook = VolRule();

      auto get_skin = [&]() {
        Range ents;
        CHKERR m_field.get_moab().get_entities_by_dimension(0, 3, ents, true);
        Skinner skin(&m_field.get_moab());
        Range skin_ents;
        CHKERR skin.find_skin(0, ents, false, skin_ents);
        // filter not owned entities, those are not on boundary
        ParallelComm *pcomm =
            ParallelComm::get_pcomm(&m_field.get_moab(), MYPCOMM_INDEX);
        Range proc_skin;
        CHKERR pcomm->filter_pstatus(skin_ents,
                                     PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                     PSTATUS_NOT, -1, &proc_skin);
        Range adj;
        for (auto d : {0, 1, 2})
          CHKERR m_field.get_moab().get_adjacencies(proc_skin, d, false, adj,
                                                    moab::Interface::UNION);
        // proc_skin.merge(adj);
        return proc_skin;
      };

      auto get_mark_skin_dofs = [&](Range &&skin) {
        auto problem_manager = m_field.getInterface<ProblemsManager>();
        auto marker_ptr = boost::make_shared<std::vector<bool>>();
        problem_manager->markDofs(simple_interface->getProblemName(), ROW, skin,
                                  *marker_ptr);
        return marker_ptr;
      };

      SmartPetscObj<Vec> v;
      CHKERR DMCreateGlobalVector_MoFEM(dm, v);
      auto mark_dofs = get_mark_skin_dofs(get_skin());

      // set operator to the volume elements
      domain_fe->getOpPtrVector().push_back(new OpVolumeSet("FIELD"));
      domain_fe->getOpPtrVector().push_back(new OpVolumeTest("FIELD"));

      domain_fe->getOpPtrVector().push_back(
          new OpSetBc("FIELD", true, mark_dofs));
      domain_fe->getOpPtrVector().push_back(new OpVolumeAssemble("FIELD", v));
      domain_fe->getOpPtrVector().push_back(new OpUnSetBc("FIELD"));

      // make integration in volume (here real calculations starts)
      CHKERR DMoFEMLoopFiniteElements(dm, simple_interface->getDomainFEName(),
                                      domain_fe);
      OpVolumeSet::entsIndices.clear();

      MOFEM_LOG_SYNCHRONISE(m_field.get_comm());
      CHKERR VecAssemblyBegin(v);
      CHKERR VecAssemblyEnd(v);
      CHKERR VecGhostUpdateBegin(v, ADD_VALUES, SCATTER_REVERSE);
      CHKERR VecGhostUpdateEnd(v, ADD_VALUES, SCATTER_REVERSE);
      CHKERR VecGhostUpdateBegin(v, INSERT_VALUES, SCATTER_FORWARD);
      CHKERR VecGhostUpdateEnd(v, INSERT_VALUES, SCATTER_FORWARD);

      CHKERR DMoFEMMeshToLocalVector(simple_interface->getDM(), v,
                                     INSERT_VALUES, SCATTER_REVERSE);

      // CHKERR VecView(v,PETSC_VIEWER_STDOUT_WORLD);

      const MoFEM::Problem *problem_ptr;
      CHKERR DMMoFEMGetProblemPtr(dm, &problem_ptr);

      double *array;
      CHKERR VecGetArray(v, &array);
      for (auto dof : *(problem_ptr->getNumeredRowDofs())) {
        auto loc_idx = dof->getPetscLocalDofIdx();

        if (dof->getFieldData() != array[loc_idx])
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "Data inconsistency");

        if ((*mark_dofs)[loc_idx]) {
          if (array[loc_idx] != 0)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                    "Dof assembled, but it should be not");
        } else {
          if (array[loc_idx] == 0)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                    "Dof not assembled, but it should be");
        }
      }

      CHKERR VecRestoreArray(v, &array);
    }
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}
