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

// namespace bio = boost::iostreams;
// using bio::tee_device;
// using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscBool flg = PETSC_TRUE;
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
    //   2.,0,0,
    //   0,2.,0,
    //   0,0,2.
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

    // Select base
    enum bases { AINSWORTH, DEMKOWICZ, LASBASETOP };
    const char *list_bases[] = {"ainsworth", "demkowicz"};

    PetscInt choice_base_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                                LASBASETOP, &choice_base_value, &flg);

    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "base not set");
    }

    FieldApproximationBase base = AINSWORTH_LEGENDRE_BASE;
    if (choice_base_value == AINSWORTH) {
      base = AINSWORTH_LEGENDRE_BASE;
    } else if (choice_base_value == DEMKOWICZ) {
      base = DEMKOWICZ_JACOBI_BASE;
    }

    // create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;
    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();

    // set entities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    CHKERR moab.create_meshset(MESHSET_SET, meshset_level0);
    
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    // fields
    CHKERR m_field.add_field("HCURL", HCURL, base, 1);
    
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "HCURL");
    
    // set app. order
    int order = 4;
    CHKERR m_field.set_field_order(root_set, MBTET, "HCURL", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "HCURL", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "HCURL", order);
    
    // build field
    CHKERR m_field.build_fields();

    // finite elements
    CHKERR m_field.add_finite_element("TET_FE");
    CHKERR m_field.add_finite_element("SKIN_FE");
    

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TET_FE", "HCURL");
    CHKERR m_field.modify_finite_element_add_field_col("TET_FE", "HCURL");
    CHKERR m_field.modify_finite_element_add_field_data("TET_FE", "HCURL");
    CHKERR m_field.modify_finite_element_add_field_row("SKIN_FE", "HCURL");
    CHKERR m_field.modify_finite_element_add_field_col("SKIN_FE", "HCURL");
    CHKERR m_field.modify_finite_element_add_field_data("SKIN_FE", "HCURL");
    
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

    Vec v;
    CHKERR m_field.getInterface<VecManager>()->vecCreateGhost("TEST_PROBLEM",
                                                              ROW, &v);
    CHKERR VecSetRandom(v, PETSC_NULL);
    CHKERR m_field.getInterface<VecManager>()->setLocalGhostVector(
        "TEST_PROBLEM", ROW, v, INSERT_VALUES, SCATTER_REVERSE);
    CHKERR VecDestroy(&v);

    struct OpTetCurl
        : public VolumeElementForcesAndSourcesCore::UserDataOperator {

      FTensor::Tensor1<double *, 3> &cUrl;
      OpTetCurl(FTensor::Tensor1<double *, 3> &curl)
          : VolumeElementForcesAndSourcesCore::UserDataOperator(
                "HCURL", UserDataOperator::OPROW),
            cUrl(curl) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBegin;

        if (data.getFieldData().size() == 0)
          MoFEMFunctionReturnHot(0);

        const unsigned int nb_gauss_pts = data.getDiffHcurlN().size1();
        const unsigned int nb_dofs = data.getFieldData().size();

        MatrixDouble curl_mat;
        FTensor::Index<'i', 3> i;

        unsigned int gg = 0;
        for (; gg < nb_gauss_pts; gg++) {
          double w = getGaussPts()(3, gg) * getVolume();
          if (getHoGaussPtsDetJac().size() == nb_gauss_pts) {
            // if ho geometry is given
            w *= getHoGaussPtsDetJac()(gg);
          }
          CHKERR getCurlOfHCurlBaseFunctions(side, type, data, gg, curl_mat);
          FTensor::Tensor1<double *, 3> t_curl(&curl_mat(0, 0), &curl_mat(0, 1),
                                               &curl_mat(0, 2), 3);
          for (unsigned int dd = 0; dd != nb_dofs; dd++) {
            cUrl(i) += w * t_curl(i) * data.getFieldData()[dd];
            ++t_curl;
          }
        }

        MoFEMFunctionReturn(0);
      }
    };

    struct MyFE : public VolumeElementForcesAndSourcesCore {

      MyFE(MoFEM::Interface &m_field)
          : VolumeElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return 2 * order+1; }; 
    };

    struct MyTriFE : public FaceElementForcesAndSourcesCore {

      MyTriFE(MoFEM::Interface &m_field)
          : FaceElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return 2 * order+1; }; 
    };

    struct OpFacesRot
        : public FaceElementForcesAndSourcesCore::UserDataOperator {

      FTensor::Tensor1<double *, 3> &cUrl;
      OpFacesRot(FTensor::Tensor1<double *, 3> &curl)
          : FaceElementForcesAndSourcesCore::UserDataOperator(
                "HCURL", UserDataOperator::OPROW),
            cUrl(curl) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBegin;

        int nb_dofs = data.getFieldData().size();
        if (nb_dofs == 0)
          MoFEMFunctionReturnHot(0);
        int nb_gauss_pts = data.getHcurlN().size1();

        auto t_curl_base = data.getFTensor1HcurlN<3>();
        // double area = getArea();
        double n0 = getNormal()[0] * 0.5;
        double n1 = getNormal()[1] * 0.5;
        double n2 = getNormal()[2] * 0.5;

        FTensor::Index<'i', 3> i;
        FTensor::Index<'j', 3> j;

        for (int gg = 0; gg < nb_gauss_pts; gg++) {
          for (int dd = 0; dd < nb_dofs; dd++) {
            double w = getGaussPts()(2, gg);
            if (getNormalsAtGaussPt().size1() == (unsigned int)nb_gauss_pts) {
              n0 = getNormalsAtGaussPt(gg)[0] * 0.5;
              n1 = getNormalsAtGaussPt(gg)[1] * 0.5;
              n2 = getNormalsAtGaussPt(gg)[2] * 0.5;
            }
            double v = data.getFieldData()[dd];
            cUrl(0) += (n1 * t_curl_base(2) - n2 * t_curl_base(1)) * w * v;
            cUrl(1) += (n2 * t_curl_base(0) - n0 * t_curl_base(2)) * w * v;
            cUrl(2) += (n0 * t_curl_base(1) - n1 * t_curl_base(0)) * w * v;
            ++t_curl_base;
          }
        }

        MoFEMFunctionReturn(0);
      }
    };

    VectorDouble curl_vol(3);
    VectorDouble curl_skin(3);
    FTensor::Tensor1<double *, 3> t_curl_vol(&curl_vol[0], &curl_vol[1],
                                             &curl_vol[2]);
    FTensor::Tensor1<double *, 3> t_curl_skin(&curl_skin[0], &curl_skin[1],
                                              &curl_skin[2]);
    curl_vol.clear();
    curl_skin.clear();

    MyFE tet_fe(m_field);
    tet_fe.getOpPtrVector().push_back(new OpTetCurl(t_curl_vol));
    MyTriFE skin_fe(m_field);
    skin_fe.getOpPtrVector().push_back(new OpFacesRot(t_curl_skin));

    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TET_FE", tet_fe);
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "SKIN_FE", skin_fe);

    std::cout.precision(12);

    std::cout << "curl_vol " << curl_vol << std::endl;
    std::cout << "curl_skin " << curl_skin << std::endl;

    FTensor::Index<'i', 3> i;
    t_curl_vol(i) -= t_curl_skin(i);
    double nrm2 = sqrt(t_curl_vol(i) * t_curl_vol(i));

    const double eps = 1e-8;
    if (fabs(nrm2) > eps) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "Curl operator not passed test\n");
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

    // mesh partitioning
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    t_curl_vol(i) = 0;
    t_curl_skin(i) = 0;

    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TET_FE", tet_fe);
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "SKIN_FE", skin_fe);

    std::cout << "curl_vol " << curl_vol << std::endl;
    std::cout << "curl_skin " << curl_skin << std::endl;

    t_curl_vol(i) -= t_curl_skin(i);
    nrm2 = sqrt(t_curl_vol(i) * t_curl_vol(i));
    if (fabs(nrm2) > eps) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "Curl operator not passed test\n");
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  
}
