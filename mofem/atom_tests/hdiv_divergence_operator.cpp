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

namespace bio = boost::iostreams;
using bio::stream;
using bio::tee_device;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    enum bases { HDIV_AINSWORTH, HDIV_DEMKOWICZ, LASTOP };

    const char *list[] = {"hdiv_ainsworth", "hdiv_demkowicz"};

    PetscBool flg;
    PetscInt choise_value = HDIV_AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list, LASTOP,
                                &choise_value, &flg);
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "base not set");
    }

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    char mesh_file_name[255];
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
#else
    CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-my_file",
                                 mesh_file_name, 255, &flg);
#endif
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    }
    const char *option;
    option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

    // //create one tet
    // double tet_coords[] = {
    //   0,0,0,
    //   0.5,0,0,
    //   0,0.5,0,
    //   0,0,0.5
    // };
    //
    // EntityHandle nodes[4];
    // for(int nn = 0;nn<4;nn++) {
    //   CHKERR moab.create_vertex(&tet_coords[3*nn],nodes[nn]);
    // }
    //
    // EntityHandle tet;
    // CHKERR moab.create_element(MBTET,nodes,4,tet);
    //
    // ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    // if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

    // create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;
    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();

    // set entitities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    CHKERR moab.create_meshset(MESHSET_SET, meshset_level0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    // fields
    switch (choise_value) {
    case HDIV_AINSWORTH:
      CHKERR m_field.add_field("HDIV", HDIV, AINSWORTH_LEGENDRE_BASE, 1);
      break;
    case HDIV_DEMKOWICZ:
      CHKERR m_field.add_field("HDIV", HDIV, DEMKOWICZ_JACOBI_BASE, 1);
      break;
    }

    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "HDIV");
    // set app. order
    int order = 5;
    CHKERR m_field.set_field_order(root_set, MBTET, "HDIV", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "HDIV", order);
    // build field
    CHKERR m_field.build_fields();

    // finite elements
    CHKERR m_field.add_finite_element("TET_FE");
    CHKERR m_field.add_finite_element("SKIN_FE");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TET_FE", "HDIV");
    CHKERR m_field.modify_finite_element_add_field_col("TET_FE", "HDIV");
    CHKERR m_field.modify_finite_element_add_field_data("TET_FE", "HDIV");
    CHKERR m_field.modify_finite_element_add_field_row("SKIN_FE", "HDIV");
    CHKERR m_field.modify_finite_element_add_field_col("SKIN_FE", "HDIV");
    CHKERR m_field.modify_finite_element_add_field_data("SKIN_FE", "HDIV");
    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET,
                                                      "TET_FE");
    Range tets;
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        BitRefLevel().set(0), BitRefLevel().set(), MBTET, tets);
    Skinner skin(&moab);
    Range skin_faces; // skin faces from 3d ents
    CHKERR skin.find_skin(0, tets, false, skin_faces);
    CHKERR m_field.add_ents_to_finite_element_by_type(skin_faces, MBTRI,
                                                      "SKIN_FE");

    // build finite elemnts
    CHKERR m_field.build_finite_elements();

    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);

    // problem
    CHKERR m_field.add_problem("TEST_PROBLEM");
    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "TET_FE");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "SKIN_FE");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);
    // build problem
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);

    // mesh partitioning
    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    struct OpTetDivergence
        : public VolumeElementForcesAndSourcesCore::UserDataOperator {

      double &dIv;
      OpTetDivergence(double &div)
          : VolumeElementForcesAndSourcesCore::UserDataOperator(
                "HDIV", UserDataOperator::OPROW),
            dIv(div) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBeginHot;

        if (type != MBTRI && type != MBTET)
          MoFEMFunctionReturnHot(0);

        if (data.getFieldData().size() == 0)
          MoFEMFunctionReturnHot(0);

        // cout << "type " << type << " side " << side << std::endl;

        int nb_gauss_pts = data.getDiffN().size1();
        int nb_dofs = data.getFieldData().size();

        VectorDouble div_vec;
        div_vec.resize(nb_dofs, 0);

        int gg = 0;
        for (; gg < nb_gauss_pts; gg++) {
          CHKERR getDivergenceOfHDivBaseFunctions(side, type, data, gg,
                                                  div_vec);
          // cout << std::fixed << div_vec << std::endl;
          unsigned int dd = 0;
          for (; dd < div_vec.size(); dd++) {
            double w = getGaussPts()(3, gg) * getVolume();
            if (getHoGaussPtsDetJac().size() > 0) {
              w *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
            }
            dIv += div_vec[dd] * w;
          }
        }

        MoFEMFunctionReturnHot(0);
      }
    };

    struct MyFE : public VolumeElementForcesAndSourcesCore {

      MyFE(MoFEM::Interface &m_field)
          : VolumeElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return order + 1; };
    };

    struct MyTriFE : public FaceElementForcesAndSourcesCore {

      MyTriFE(MoFEM::Interface &m_field)
          : FaceElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return order + 1; };
    };

    struct OpFacesFluxes
        : public FaceElementForcesAndSourcesCore::UserDataOperator {

      double &dIv;
      OpFacesFluxes(double &div)
          : FaceElementForcesAndSourcesCore::UserDataOperator(
                "HDIV", UserDataOperator::OPROW),
            dIv(div) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBeginHot;

        if (type != MBTRI)
          MoFEMFunctionReturnHot(0);

        int nb_gauss_pts = data.getN().size1();
        int nb_dofs = data.getFieldData().size();

        int gg = 0;
        for (; gg < nb_gauss_pts; gg++) {
          int dd = 0;
          for (; dd < nb_dofs; dd++) {
            double area;
            VectorDouble n;
            if (getNormalsAtGaussPts().size1() == (unsigned int)nb_gauss_pts) {
              n = getNormalsAtGaussPts(gg);
              area = norm_2(getNormalsAtGaussPts(gg)) * 0.5;
            } else {
              n = getNormal();
              area = getArea();
            }
            n /= norm_2(n);
            dIv += (n[0] * data.getVectorN<3>(gg)(dd, 0) +
                    n[1] * data.getVectorN<3>(gg)(dd, 1) +
                    n[2] * data.getVectorN<3>(gg)(dd, 2)) *
                   getGaussPts()(2, gg) * area;
          }
          // cout << getNormal() << std::endl;
        }

        MoFEMFunctionReturnHot(0);
      }
    };

    double divergence_vol = 0;
    double divergence_skin = 0;

    MyFE tet_fe(m_field);
    tet_fe.getOpPtrVector().push_back(new OpTetDivergence(divergence_vol));

    MyTriFE skin_fe(m_field);
    skin_fe.getOpPtrVector().push_back(new OpFacesFluxes(divergence_skin));

    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TET_FE", tet_fe);
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "SKIN_FE", skin_fe);

    std::cout.precision(12);

    std::cout << "divergence_vol " << divergence_vol << std::endl;
    std::cout << "divergence_skin " << divergence_skin << std::endl;

    const double eps = 1e-6;
    if (fabs(divergence_skin - divergence_vol) > eps) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
               "invalid surface flux or divergence or both\n", divergence_skin,
               divergence_vol);
    }

    CHKERR m_field.add_field("MESH_NODE_POSITIONS", H1, AINSWORTH_LEGENDRE_BASE,
                             3);

    CHKERR m_field.add_ents_to_field_by_type(0, MBTET, "MESH_NODE_POSITIONS");
    CHKERR m_field.set_field_order(0, MBVERTEX, "MESH_NODE_POSITIONS", 1);
    CHKERR m_field.set_field_order(0, MBEDGE, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(0, MBTRI, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(0, MBTET, "MESH_NODE_POSITIONS", 2);

    CHKERR m_field.modify_finite_element_add_field_data("TET_FE",
                                                        "MESH_NODE_POSITIONS");
    CHKERR m_field.modify_finite_element_add_field_data("SKIN_FE",
                                                        "MESH_NODE_POSITIONS");

    CHKERR m_field.build_fields();
    // project geometry form 10 node tets on higher order approx. functions
    Projection10NodeCoordsOnField ent_method(m_field, "MESH_NODE_POSITIONS");
    CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method);

    CHKERR m_field.build_finite_elements();
    CHKERR m_field.build_adjacencies(bit_level0);

    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);

    // mesh partitioning
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    // for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(m_field,"MESH_NODE_POSITIONS",MBVERTEX,dof))
    // {
    //   EntityHandle vert = (*dof)->getEnt();
    //   double coords[3];
    //   CHKERR moab.get_coords(&vert,1,coords);
    //   coords[0] *= 2;
    //   coords[1] *= 4;
    //   coords[2] *= 0.5;
    //
    //   (*dof)->getFieldData() = coords[(*dof)->getDofCoeffIdx()];
    // }

    divergence_vol = 0;
    divergence_skin = 0;

    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TET_FE", tet_fe);
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "SKIN_FE", skin_fe);

    std::cout.precision(12);

    std::cout << "divergence_vol " << divergence_vol << std::endl;
    std::cout << "divergence_skin " << divergence_skin << std::endl;

    if (fabs(divergence_skin - divergence_vol) > eps) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
               "invalid surface flux or divergence or both\n", divergence_skin,
               divergence_vol);
    }
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
