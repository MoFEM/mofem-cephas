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

static const double face_coords[4][9] = {{0, 0, 0, 1, 0, 0, 0, 0, 1},
                                         {1, 0, 0, 0, 1, 0, 0, 0, 1},
                                         {0, 0, 0, 0, 1, 0, 0, 0, 1},
                                         {0, 0, 0, 1, 0, 0, 0, 1, 0}};

static const double edge_coords[6][6] = {
    {0, 0, 0, 1, 0, 0}, {1, 0, 0, 0, 1, 0}, {0, 1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 1}, {1, 0, 0, 0, 0, 1}, {0, 1, 0, 0, 0, 1}};

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
    ierr = PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
    CHKERRG(ierr);
#else
    ierr = PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-my_file",
                                 mesh_file_name, 255, &flg);
    CHKERRG(ierr);
#endif
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    }
    const char *option;
    option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // set entitities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    CHKERR moab.create_meshset(MESHSET_SET, meshset_level0);
    
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);
    

    // Fields
    CHKERR m_field.add_field("MESH_NODE_POSITIONS", H1, AINSWORTH_LEGENDRE_BASE,
                             3);
    
    CHKERR m_field.add_field("H1", H1, AINSWORTH_LEGENDRE_BASE, 1);
    

    // FE
    CHKERR m_field.add_finite_element("TET_FE");
    CHKERR m_field.add_finite_element("TRI_FE");
    CHKERR m_field.add_finite_element("SKIN_FE");
    CHKERR m_field.add_finite_element("EDGE_FE");
    

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TET_FE", "H1");
    CHKERR m_field.modify_finite_element_add_field_col("TET_FE", "H1");
    CHKERR m_field.modify_finite_element_add_field_data("TET_FE", "H1");
    CHKERR m_field.modify_finite_element_add_field_data("TET_FE",
                                                        "MESH_NODE_POSITIONS");
    

    CHKERR m_field.modify_finite_element_add_field_row("SKIN_FE", "H1");
    CHKERR m_field.modify_finite_element_add_field_col("SKIN_FE", "H1");
    CHKERR m_field.modify_finite_element_add_field_data("SKIN_FE", "H1");
    CHKERR m_field.modify_finite_element_add_field_data("SKIN_FE",
                                                        "MESH_NODE_POSITIONS");
    

    CHKERR m_field.modify_finite_element_add_field_row("TRI_FE", "H1");
    CHKERR m_field.modify_finite_element_add_field_col("TRI_FE", "H1");
    CHKERR m_field.modify_finite_element_add_field_data("TRI_FE", "H1");
    CHKERR m_field.modify_finite_element_add_field_data("TRI_FE",
                                                        "MESH_NODE_POSITIONS");
    

    CHKERR m_field.modify_finite_element_add_field_row("EDGE_FE", "H1");
    CHKERR m_field.modify_finite_element_add_field_col("EDGE_FE", "H1");
    CHKERR m_field.modify_finite_element_add_field_data("EDGE_FE", "H1");
    CHKERR m_field.modify_finite_element_add_field_data("EDGE_FE",
                                                        "MESH_NODE_POSITIONS");

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");
    

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "TET_FE");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "SKIN_FE");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "TRI_FE");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "EDGE_FE");
    
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);
    

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "H1");
    

    // add entities to finite element
    CHKERR
        m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "TET_FE");
    

    Range tets;
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        BitRefLevel().set(0), BitRefLevel().set(), MBTET, tets);
    
    Skinner skin(&moab);
    Range skin_faces; // skin faces from 3d ents
    CHKERR skin.find_skin(0, tets, false, skin_faces);
    CHKERR m_field.add_ents_to_finite_element_by_type(skin_faces, MBTRI,
                                                      "SKIN_FE");

    Range faces;
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        BitRefLevel().set(0), BitRefLevel().set(), MBTRI, faces);
    
    faces = subtract(faces, skin_faces);
    CHKERR m_field.add_ents_to_finite_element_by_type(faces, MBTRI, "TRI_FE");

    Range edges;
    CHKERR moab.get_adjacencies(faces, 1, false, edges, moab::Interface::UNION);
    CHKERR m_field.add_ents_to_finite_element_by_type(edges, MBEDGE, "EDGE_FE");
    
    Range skin_edges;
    CHKERR moab.get_adjacencies(skin_faces, 1, false, skin_edges,
                                moab::Interface::UNION);
    CHKERR m_field.add_ents_to_finite_element_by_type(edges, MBEDGE, "EDGE_FE");
    

    // set app. order
    int order = 4;
    CHKERR m_field.set_field_order(root_set, MBTET, "H1", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "H1", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "H1", order);

    CHKERR m_field.add_ents_to_field_by_type(0, MBTET, "MESH_NODE_POSITIONS");
    CHKERR m_field.set_field_order(0, MBVERTEX, "MESH_NODE_POSITIONS", 1);
    CHKERR m_field.set_field_order(0, MBEDGE, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(0, MBTRI, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(0, MBTET, "MESH_NODE_POSITIONS", 2);
    
    /****/
    // build database
    // build field
    CHKERR m_field.build_fields();
    
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);
    

    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    
    // build problem
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);
    

    // project geometry form 10 node tets on higher order approx. functions
    Projection10NodeCoordsOnField ent_method(m_field, "MESH_NODE_POSITIONS");
    CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method);

    /****/
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
    

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;
    std::ofstream ofs("forces_and_sources_hdiv_continuity_check.txt");
    TeeDevice my_tee(std::cout, ofs);
    TeeStream my_split(my_tee);

    struct OpTetFluxes
        : public VolumeElementForcesAndSourcesCore::UserDataOperator {

      MoFEM::Interface &mField;
      Tag tH;

      OpTetFluxes(MoFEM::Interface &m_field, Tag th)
          : VolumeElementForcesAndSourcesCore::UserDataOperator(
                "H1", UserDataOperator::OPROW),
            mField(m_field), tH(th) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
        MoFEMFunctionBegin;

        if (data.getFieldData().size() == 0)
          MoFEMFunctionReturnHot(0);

        if (type == MBTRI) {

          boost::shared_ptr<const NumeredEntFiniteElement> mofem_fe =
              getNumeredEntFiniteElementPtr();
          SideNumber_multiIndex &side_table = mofem_fe->getSideNumberTable();
          EntityHandle face = side_table.get<1>()
                                  .find(boost::make_tuple(type, side))
                                  ->get()
                                  ->ent;
          int sense = side_table.get<1>()
                          .find(boost::make_tuple(type, side))
                          ->get()
                          ->sense;

          double t = 0;
          int nb_dofs = data.getN().size2();
          for (int dd = 0; dd < nb_dofs; dd++) {
            t += data.getN(side)(dd) * data.getFieldData()[dd];
          }

          double *t_ptr;
          CHKERR mField.get_moab().tag_get_by_ptr(tH, &face, 1,
                                                  (const void **)&t_ptr);
          
          *t_ptr += sense * t;
        }

        if (type == MBEDGE) {

          boost::shared_ptr<const NumeredEntFiniteElement> mofem_fe =
              getNumeredEntFiniteElementPtr();
          SideNumber_multiIndex &side_table = mofem_fe->getSideNumberTable();
          EntityHandle edge = side_table.get<1>()
                                  .find(boost::make_tuple(type, side))
                                  ->get()
                                  ->ent;
          Range adj_tets;
          CHKERR mField.get_moab().get_adjacencies(&edge, 1, 3, false, adj_tets,
                                                   moab::Interface::UNION);
          const int nb_adj_tets = adj_tets.size();

          double t = 0;
          int nb_dofs = data.getN().size2();
          for (int dd = 0; dd < nb_dofs; dd++) {
            t += data.getN(4 + side)(dd) * data.getFieldData()[dd];
          }

          double *t_ptr;
          CHKERR mField.get_moab().tag_get_by_ptr(tH, &edge, 1,
                                                  (const void **)&t_ptr);
          
          *t_ptr += t / nb_adj_tets;
        }

        MoFEMFunctionReturn(0);
      }
    };

    struct MyTetFE : public VolumeElementForcesAndSourcesCore {

      MyTetFE(MoFEM::Interface &m_field)
          : VolumeElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return -1; };

      MatrixDouble N_tri;
      MoFEMErrorCode setGaussPts(int order) {

        MoFEMFunctionBegin;

        try {

          N_tri.resize(1, 3);
          CHKERR ShapeMBTRI(&N_tri(0, 0), G_TRI_X1, G_TRI_Y1, 1);
          

          gaussPts.resize(4, 4 + 6);
          int ff = 0;
          for (; ff < 4; ff++) {
            int dd = 0;
            for (; dd < 3; dd++) {
              gaussPts(dd, ff) =
                  cblas_ddot(3, &N_tri(0, 0), 1, &face_coords[ff][dd], 3);
            }
            gaussPts(3, ff) = G_TRI_W1[0];
          }

          int ee = 0;
          for (; ee < 6; ee++) {
            int dd = 0;
            for (; dd < 3; dd++) {
              gaussPts(dd, 4 + ee) =
                  (edge_coords[ee][0 + dd] + edge_coords[ee][3 + dd]) * 0.5;
            }
            gaussPts(3, 4 + ee) = 1;
          }

          // std::cerr << gaussPts << std::endl;

        } catch (std::exception &ex) {
          std::ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__
             << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF, MOFEM_STD_EXCEPTION_THROW, ss.str().c_str());
        }

        MoFEMFunctionReturn(0);
      }
    };

    struct OpFacesSkinFluxes
        : public FaceElementForcesAndSourcesCore::UserDataOperator {

      MoFEM::Interface &mField;
      Tag tH1, tH2;
      TeeStream &mySplit;

      OpFacesSkinFluxes(MoFEM::Interface &m_field, Tag th1, Tag th2,
                        TeeStream &my_split)
          : FaceElementForcesAndSourcesCore::UserDataOperator(
                "H1", UserDataOperator::OPROW),
            mField(m_field), tH1(th1), tH2(th2), mySplit(my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
        MoFEMFunctionBegin;

        if (type != MBTRI)
          MoFEMFunctionReturnHot(0);
        EntityHandle face = getNumeredEntFiniteElementPtr()->getEnt();

        double *t_ptr;
        CHKERR mField.get_moab().tag_get_by_ptr(tH1, &face, 1,
                                                (const void **)&t_ptr);
        
        double *tn_ptr;
        CHKERR mField.get_moab().tag_get_by_ptr(tH2, &face, 1,
                                                (const void **)&tn_ptr);
        

        *tn_ptr = *t_ptr;

        int dd = 0;
        int nb_dofs = data.getN().size2();
        for (; dd < nb_dofs; dd++) {
          double val = data.getFieldData()[dd];
          *tn_ptr += -data.getN()(0, dd) * val;
        }

        const double eps = 1e-8;
        if (fabs(*tn_ptr) > eps) {
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "H1 continuity failed %6.4e", *tn_ptr);
        }

        mySplit.precision(5);
        mySplit << face << " " << /*std::fixed <<*/ fabs(*tn_ptr) << std::endl;

        MoFEMFunctionReturn(0);
      }
    };

    struct OpFacesFluxes
        : public FaceElementForcesAndSourcesCore::UserDataOperator {

      MoFEM::Interface &mField;
      Tag tH1, tH2;
      TeeStream &mySplit;

      OpFacesFluxes(MoFEM::Interface &m_field, Tag _th1, Tag _th2,
                    TeeStream &my_split)
          : FaceElementForcesAndSourcesCore::UserDataOperator(
                "H1", UserDataOperator::OPROW),
            mField(m_field), tH1(_th1), tH2(_th2), mySplit(my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
        MoFEMFunctionBegin;

        if (type != MBTRI)
          MoFEMFunctionReturnHot(0);

        EntityHandle face = getNumeredEntFiniteElementPtr()->getEnt();

        double *t_ptr;
        CHKERR mField.get_moab().tag_get_by_ptr(tH1, &face, 1,
                                                (const void **)&t_ptr);
        
        double *tn_ptr;
        CHKERR mField.get_moab().tag_get_by_ptr(tH2, &face, 1,
                                                (const void **)&tn_ptr);
        

        *tn_ptr = *t_ptr;

        const double eps = 1e-8;
        if (fabs(*tn_ptr) > eps) {
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "H1 continuity failed %6.4e", *tn_ptr);
        }

        mySplit.precision(5);

        mySplit << face << " " << /*std::fixed <<*/ fabs(*tn_ptr) << std::endl;

        MoFEMFunctionReturn(0);
      }
    };

    struct MyTriFE : public FaceElementForcesAndSourcesCore {

      MyTriFE(MoFEM::Interface &m_field)
          : FaceElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return -1; };

      MoFEMErrorCode setGaussPts(int order) {
        MoFEMFunctionBegin;

        gaussPts.resize(3, 1);
        gaussPts(0, 0) = G_TRI_X1[0];
        gaussPts(1, 0) = G_TRI_Y1[0];
        gaussPts(2, 0) = G_TRI_W1[0];

        MoFEMFunctionReturn(0);
      }
    };

    struct OpEdgesFluxes
        : public EdgeElementForcesAndSourcesCore::UserDataOperator {

      MoFEM::Interface &mField;
      Tag tH1, tH2;
      TeeStream &mySplit;

      OpEdgesFluxes(MoFEM::Interface &m_field, Tag _th1, Tag _th2,
                    TeeStream &my_split)
          : EdgeElementForcesAndSourcesCore::UserDataOperator(
                "H1", UserDataOperator::OPROW),
            mField(m_field), tH1(_th1), tH2(_th2), mySplit(my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
        MoFEMFunctionBegin;

        if (type != MBEDGE)
          MoFEMFunctionReturnHot(0);

        EntityHandle edge = getNumeredEntFiniteElementPtr()->getEnt();

        double *t_ptr;
        CHKERR mField.get_moab().tag_get_by_ptr(tH1, &edge, 1,
                                                (const void **)&t_ptr);
        
        double *tn_ptr;
        CHKERR mField.get_moab().tag_get_by_ptr(tH2, &edge, 1,
                                                (const void **)&tn_ptr);
        

        *tn_ptr = *t_ptr;

        double tn = 0;
        int nb_dofs = data.getN().size2();
        int dd = 0;
        for (; dd < nb_dofs; dd++) {
          double val = data.getFieldData()[dd];
          tn += data.getN()(0, dd) * val;
        }

        // mySplit << *tn_ptr << "  " << tn << " " << getLength() << endl;
        *tn_ptr -= tn;

        // mySplit << getTangetAtGaussPts() << " " << getDirection() << endl;

        // cerr << t_ptr[0] << " " << t_ptr[1] << " " << t_ptr[2] << endl;

        // const double eps = 1e-8;
        // if(fabs(*tn_ptr)>eps) {
        //   SETERRQ1(
        //     PETSC_COMM_SELF,
        //     MOFEM_ATOM_TEST_INVALID,
        //     "H1 continuity failed %6.4e",
        //     *tn_ptr
        //   );
        // }

        mySplit.precision(5);

        mySplit << edge << " " << /*std::fixed <<*/ fabs(*tn_ptr) << std::endl;

        MoFEMFunctionReturn(0);
      }
    };

    struct MyEdgeFE : public EdgeElementForcesAndSourcesCore {

      MyEdgeFE(MoFEM::Interface &m_field)
          : EdgeElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return -1; };

      MoFEMErrorCode setGaussPts(int order) {
        MoFEMFunctionBegin;

        gaussPts.resize(2, 1);
        gaussPts(0, 0) = 0.5;
        gaussPts(1, 0) = 1;

        MoFEMFunctionReturn(0);
      }
    };

    MyTetFE tet_fe(m_field);
    MyTriFE tri_fe(m_field);
    MyTriFE skin_fe(m_field);
    MyEdgeFE edge_fe(m_field);

    Tag th1;
    double def_val[] = {0, 0, 0};
    CHKERR moab.tag_get_handle("T", 3, MB_TYPE_DOUBLE, th1,
                               MB_TAG_CREAT | MB_TAG_SPARSE, &def_val);

    tet_fe.getOpPtrVector().push_back(
        new OpCalculateHOCoords("MESH_NODE_POSITIONS"));
    tet_fe.getOpPtrVector().push_back(new OpTetFluxes(m_field, th1));

    Tag th2;
    CHKERR moab.tag_get_handle("TN", 1, MB_TYPE_DOUBLE, th2,
                               MB_TAG_CREAT | MB_TAG_SPARSE, &def_val);

    tri_fe.getOpPtrVector().push_back(
        new OpCalculateHOCoords("MESH_NODE_POSITIONS"));
    tri_fe.getOpPtrVector().push_back(
        new OpGetHONormalsOnFace("MESH_NODE_POSITIONS"));
    tri_fe.getOpPtrVector().push_back(
        new OpFacesFluxes(m_field, th1, th2, my_split));
    skin_fe.getOpPtrVector().push_back(
        new OpFacesSkinFluxes(m_field, th1, th2, my_split));

    edge_fe.getOpPtrVector().push_back(
        new OpGetHOTangentsOnEdge("MESH_NODE_POSITIONS"));
    edge_fe.getOpPtrVector().push_back(
        new OpCalculateHOCoords("MESH_NODE_POSITIONS"));
    edge_fe.getOpPtrVector().push_back(
        new OpEdgesFluxes(m_field, th1, th2, my_split));

    for (Range::iterator fit = faces.begin(); fit != faces.end(); fit++) {
      CHKERR moab.tag_set_data(th1, &*fit, 1, &def_val);
      
      CHKERR moab.tag_set_data(th2, &*fit, 1, &def_val);
      
    }

    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TET_FE", tet_fe);
    my_split << "internal\n";
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TRI_FE", tri_fe);
    my_split << "skin\n";
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "SKIN_FE", skin_fe);
    my_split << "edges\n";
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "EDGE_FE", edge_fe);
    

    EntityHandle meshset;
    CHKERR moab.create_meshset(MESHSET_SET, meshset);
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        BitRefLevel().set(0), BitRefLevel().set(), MBTRI, meshset);
    
    CHKERR moab.write_file("out.vtk", "VTK", "", &meshset, 1);
    
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}
