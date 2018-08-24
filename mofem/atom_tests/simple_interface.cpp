/**
 * \file simple_interface.cpp
 * \ingroup mofem_simple_interface
 * \example simple_interface.cpp
 *
 * Calculate volume by integrating volume elements and using divergence theorem
 * by integrating surface elements.
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

struct OpVolume : public VolumeElementForcesAndSourcesCore::UserDataOperator {
  Vec vOl;
  OpVolume(const std::string &field_name, Vec vol)
      : VolumeElementForcesAndSourcesCore::UserDataOperator(field_name, OPROW),
        vOl(vol) {}
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {

    MoFEMFunctionBegin;
    if (type != MBVERTEX)
      MoFEMFunctionReturnHot(0);
    const int nb_int_pts = getGaussPts().size2();
    // cerr << nb_int_pts << endl;
    auto t_w = getFTensor0IntegrationWeight();
    auto t_ho_det = getFTenosr0HoMeasure();
    double v = getMeasure();
    double vol = 0;
    for (int gg = 0; gg != nb_int_pts; gg++) {
      vol += t_w * t_ho_det * v;
      // cerr << t_ho_det << endl;
      ++t_w;
      ++t_ho_det;
    }
    CHKERR VecSetValue(vOl, 0, vol, ADD_VALUES);
    MoFEMFunctionReturn(0);
  }
  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        DataForcesAndSourcesCore::EntData &row_data,
                        DataForcesAndSourcesCore::EntData &col_data) {
    MoFEMFunctionBeginHot;
    // PetscPrintf(PETSC_COMM_WORLD,"domain: calculate matrix\n");
    MoFEMFunctionReturnHot(0);
  }
};

struct OpFace : public FaceElementForcesAndSourcesCore::UserDataOperator {
  Vec vOl;
  OpFace(const std::string &field_name, Vec vol)
      : FaceElementForcesAndSourcesCore::UserDataOperator(field_name, OPROW),
        vOl(vol) {}
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {

    MoFEMFunctionBegin;
    if (type != MBVERTEX)
      MoFEMFunctionReturnHot(0);
    const int nb_int_pts = getGaussPts().size2();
    auto t_normal = getFTensor1NormalsAtGaussPts();
    auto t_w = getFTensor0IntegrationWeight();
    auto t_coords = getFTensor1HoCoordsAtGaussPts();
    FTensor::Index<'i', 3> i;
    double vol = 0;
    for (int gg = 0; gg != nb_int_pts; gg++) {
      vol += (t_coords(i) * t_normal(i)) * t_w;
      ++t_normal;
      ++t_w;
      ++t_coords;
    }
    vol /= 6;
    CHKERR VecSetValue(vOl, 0, vol, ADD_VALUES);
    MoFEMFunctionReturn(0);
  }
  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        DataForcesAndSourcesCore::EntData &row_data,
                        DataForcesAndSourcesCore::EntData &col_data) {
    MoFEMFunctionBeginHot;
    // PetscPrintf(PETSC_COMM_WORLD,"boundary: calculate matrix\n");
    MoFEMFunctionReturnHot(0);
  }
};

struct VolRule {
  int operator()(int, int, int) const { return 2; }
};
struct FaceRule {
  int operator()(int, int, int) const { return 4; }
};

int main(int argc, char *argv[]) {

  //

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
      CHKERR simple_interface->addDomainField("MESH_NODE_POSITIONS", H1,
                                              AINSWORTH_LEGENDRE_BASE, 3);
      CHKERR simple_interface->addBoundaryField("MESH_NODE_POSITIONS", H1,
                                                AINSWORTH_LEGENDRE_BASE, 3);
      // set fields order
      CHKERR simple_interface->setFieldOrder("MESH_NODE_POSITIONS", 1);
      // setup problem
      CHKERR simple_interface->setUp();
      // Project mesh coordinate on mesh
      Projection10NodeCoordsOnField ent_method(m_field, "MESH_NODE_POSITIONS");
      CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method);
      DM dm;
      // get dm
      CHKERR simple_interface->getDM(&dm);
      CHKERRG(ierr);
      // create elements
      boost::shared_ptr<ForcesAndSourcesCore> domain_fe =
          boost::shared_ptr<ForcesAndSourcesCore>(
              new VolumeElementForcesAndSourcesCore(m_field));
      boost::shared_ptr<ForcesAndSourcesCore> boundary_fe =
          boost::shared_ptr<ForcesAndSourcesCore>(
              new FaceElementForcesAndSourcesCore(m_field));
      // set integration rule
      domain_fe->getRuleHook = VolRule();
      boundary_fe->getRuleHook = FaceRule();
      // create distributed vector to accumulate values from processors.
      int ghosts[] = {0};
      Vec vol, surf_vol;
      CHKERR VecCreateGhost(PETSC_COMM_WORLD,
                            m_field.get_comm_rank() == 0 ? 1 : 0, 1,
                            m_field.get_comm_rank() == 0 ? 0 : 1, ghosts, &vol);
      CHKERR VecDuplicate(vol, &surf_vol);
      // set operator to the volume element
      domain_fe->getOpPtrVector().push_back(
          new OpVolume("MESH_NODE_POSITIONS", vol));
      // set operator to the face element
      boundary_fe->getOpPtrVector().push_back(
          new OpFace("MESH_NODE_POSITIONS", surf_vol));
      // make integration in volume (here real calculations starts)
      CHKERR DMoFEMLoopFiniteElements(dm, simple_interface->getDomainFEName(),
                                      domain_fe);
      // make integration on boundary
      CHKERR DMoFEMLoopFiniteElements(dm, simple_interface->getBoundaryFEName(),
                                      boundary_fe);
      // assemble volumes from processors and accumulate on processor of rank 0
      CHKERR VecAssemblyBegin(vol);
      CHKERR VecAssemblyEnd(vol);
      CHKERR VecGhostUpdateBegin(vol, ADD_VALUES, SCATTER_REVERSE);
      CHKERR VecGhostUpdateEnd(vol, ADD_VALUES, SCATTER_REVERSE);
      CHKERR VecAssemblyBegin(surf_vol);
      CHKERR VecAssemblyEnd(surf_vol);
      CHKERR VecGhostUpdateBegin(surf_vol, ADD_VALUES, SCATTER_REVERSE);
      CHKERR VecGhostUpdateEnd(surf_vol, ADD_VALUES, SCATTER_REVERSE);
      if (m_field.get_comm_rank() == 0) {
        double *a_vol;
        CHKERR VecGetArray(vol, &a_vol);
        double *a_surf_vol;
        CHKERR VecGetArray(surf_vol, &a_surf_vol);
        cout << "Volume = " << a_vol[0] << endl;
        cout << "Surf Volume = " << a_surf_vol[0] << endl;
        if (fabs(a_vol[0] - a_surf_vol[0]) > 1e-12) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Should be zero");
        }
        CHKERR VecRestoreArray(vol, &a_vol);
        CHKERR VecRestoreArray(vol, &a_surf_vol);
      }
      // destroy vector
      CHKERR VecDestroy(&vol);
      CHKERR VecDestroy(&surf_vol);
      // destroy dm
      CHKERR DMDestroy(&dm);
    }
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}
