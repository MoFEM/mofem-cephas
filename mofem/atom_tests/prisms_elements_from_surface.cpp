/** \file prisms_elements_from_surface.cpp
  \example prisms_elements_from_surface.cpp
  \brief Adding prims on the surface

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

namespace bio = boost::iostreams;
using bio::stream;
using bio::tee_device;

using namespace MoFEM;

static char help[] = "...\n\n";
static int debug = 1;

static constexpr int precision_exponent = 4;
static constexpr int number_of_prisms_layers = 18;
static constexpr double delta =
    1. / static_cast<double>(number_of_prisms_layers);
static constexpr std::array<double, 3> d3 = {0, 0, 0};
static constexpr std::array<double, 3> d4 = {0, 0, delta};

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Read parameters from line command
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
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
              "error -my_file (MESH FILE NEEDED)");
    }

    // Read mesh to MOAB
    const char *option;
    option = ""; //"PARALLEL=BCAST;";//;DEBUG_IO";
    CHKERR moab.load_file(mesh_file_name, 0, option);
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

    Range tris;
    CHKERR moab.get_entities_by_type(0, MBTRI, tris, false);
    Range tris_verts;
    CHKERR moab.get_connectivity(tris, tris_verts);
    MatrixDouble tri_coords(tris_verts.size(), 3);
    CHKERR moab.get_coords(tris_verts, &tri_coords(0, 0));

    // Create MoFEM (Joseph) databas
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    PrismsFromSurfaceInterface *prisms_from_surface_interface;
    CHKERR m_field.getInterface(prisms_from_surface_interface);

    Range prisms;
    CHKERR prisms_from_surface_interface->createPrisms(tris, prisms);
    prisms_from_surface_interface->setThickness(prisms, d3.data(), d4.data());
    Range add_prims_layer;
    Range extrude_prisms = prisms;

    for (int ll = 1; ll != number_of_prisms_layers; ++ll) {
      prisms_from_surface_interface->createdVertices.clear();
      CHKERR prisms_from_surface_interface->createPrismsFromPrisms(
          extrude_prisms, false, add_prims_layer);
      prisms_from_surface_interface->setThickness(add_prims_layer, d3.data(),
                                                  d4.data());
      extrude_prisms = add_prims_layer;
      prisms.merge(add_prims_layer);
      add_prims_layer.clear();
    }

    struct CoordsAndHandle {

      inline static double getArg(double x) {
        return std::round(x * pow(10., precision_exponent - 1));
      };

      int x, y, z;
      EntityHandle node;
      CoordsAndHandle(const double *coords, EntityHandle v)
          : x(getArg(coords[0])), y(getArg(coords[1])), z(getArg(coords[2])),
            node(v) {}
    };

    typedef multi_index_container<
        CoordsAndHandle,
        indexed_by<

            hashed_unique<composite_key<
                CoordsAndHandle,
                member<CoordsAndHandle, int, &CoordsAndHandle::x>,
                member<CoordsAndHandle, int, &CoordsAndHandle::y>,
                member<CoordsAndHandle, int, &CoordsAndHandle::z>>>

            >>
        MapCoords;

    EntityHandle meshset;
    CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, meshset);
    CHKERR moab.add_entities(meshset, prisms);

    MapCoords map_coords;
    Range verts;
    CHKERR moab.get_connectivity(prisms, verts);
    MatrixDouble coords(verts.size(), 3);
    CHKERR moab.get_coords(verts, &coords(0, 0));

    for (size_t v = 0; v != verts.size(); ++v)
      map_coords.insert(CoordsAndHandle(&coords(v, 0), verts[v]));

    EntityHandle one_prism_meshset;
    CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER,
                               one_prism_meshset);

    std::array<double, 18> one_prism_coords = {0, 0, 0, 1, 0, 0, 0, 1, 0,
                                               0, 0, 1, 1, 0, 1, 0, 1, 1};
    std::array<EntityHandle, 6> one_prism_nodes;
    for (int n = 0; n != 6; ++n)
      CHKERR moab.create_vertex(&one_prism_coords[3 * n], one_prism_nodes[n]);
    EntityHandle one_prism;
    CHKERR m_field.get_moab().create_element(MBPRISM, one_prism_nodes.data(), 6,
                                             one_prism);
    Range one_prism_range;
    one_prism_range.insert(one_prism);
    CHKERR moab.add_entities(one_prism_meshset, one_prism_range);
    Range one_prism_verts;
    CHKERR moab.get_connectivity(one_prism_range, one_prism_verts);
    CHKERR moab.add_entities(one_prism_meshset, one_prism_verts);
    Range one_prism_adj_ents;
    for (int d = 1; d != 3; ++d)
      CHKERR moab.get_adjacencies(one_prism_range, d, true, one_prism_adj_ents,
                                  moab::Interface::UNION);
    CHKERR moab.add_entities(one_prism_meshset, one_prism_adj_ents);

    cerr << one_prism_adj_ents << endl;

    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setEntitiesBitRefLevel(
        one_prism_range, bit_level0);
    CHKERR prisms_from_surface_interface->seedPrismsEntities(one_prism_range,
                                                             bit_level0);

    // Fields
    CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.add_ents_to_field_by_type(one_prism_meshset, MBPRISM,
                                             "FIELD1", VERY_NOISY);

    CHKERR m_field.set_field_order(one_prism_meshset, MBVERTEX, "FIELD1", 1);
    CHKERR m_field.set_field_order(one_prism_meshset, MBEDGE, "FIELD1", 3,
                                   VERY_NOISY);
    CHKERR m_field.set_field_order(one_prism_meshset, MBTRI, "FIELD1", 3);
    CHKERR m_field.set_field_order(one_prism_meshset, MBQUAD, "FIELD1", 5,
                                   VERY_NOISY);
    CHKERR m_field.set_field_order(one_prism_meshset, MBPRISM, "FIELD1", 7,
                                   VERY_NOISY);
    CHKERR m_field.build_fields(VERY_NOISY);

    // FE
    CHKERR m_field.add_finite_element("TEST_FE1");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE1", "FIELD1");

    CHKERR m_field.add_ents_to_finite_element_by_type(one_prism_range, MBPRISM,
                                                      "TEST_FE1");

    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // //build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM",
                                                     "TEST_FE1");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);

    // build problem
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);
    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;

    struct MyOp : public FatPrismElementForcesAndSourcesCore::UserDataOperator {

      MyOp(moab::Interface &post_proc, MapCoords &map_coords)
          : FatPrismElementForcesAndSourcesCore::UserDataOperator(
                "FIELD1", "FIELD1",
                ForcesAndSourcesCore::UserDataOperator::OPROW),
            postProc(post_proc), mapCoords(map_coords) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        constexpr double def_val[] = {0, 0, 0};
        MoFEMFunctionBegin;
        switch (type) {
        case MBVERTEX:
        case MBEDGE:
        case MBTRI:
        case MBQUAD:
        case MBPRISM:
          break;
        default:
          MoFEMFunctionReturnHot(0);
        }
        if (type == MBTRI && (side != 3 && side != 4))
          MoFEMFunctionReturnHot(0);
        if (type == MBQUAD && (side == 3 || side == 4))
          MoFEMFunctionReturnHot(0);

        if (type == MBVERTEX) {
          const size_t nb_gauss_pts = getGaussPts().size2();
          auto &coords_at_pts = getGaussPts();
          nodeHandles.reserve(nb_gauss_pts);
          nodeHandles.clear();
          for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
            auto t = boost::make_tuple(
                CoordsAndHandle::getArg(coords_at_pts(0, gg)),
                CoordsAndHandle::getArg(coords_at_pts(1, gg)),
                CoordsAndHandle::getArg(coords_at_pts(2, gg)));

            auto it = mapCoords.find(t);
            if (it != mapCoords.end())
              nodeHandles.emplace_back(it->node);
            else
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "Vertex not found");
          }

          MatrixDouble node_coords(nb_gauss_pts, 3);
          CHKERR postProc.get_coords(&nodeHandles[0], nodeHandles.size(),
                                     &node_coords(0, 0));
          for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
            constexpr double eps = 1e-12;
            for (auto d : {0, 1, 2})
              if (std::abs(node_coords(gg, d) - getGaussPts()(d, gg)) > eps) {
                SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                        "Inconsistency between node coords and integration "
                        "points");
              }
            for (auto d : {0, 1, 2})
              if (std::abs(node_coords(gg, d) - getCoordsAtGaussPts()(gg, d)) >
                  eps)
                SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                        "Inconsistency between node coords and integration "
                        "points");
          }

          Tag th;
          CHKERR postProc.tag_get_handle("Coords", 3, MB_TYPE_DOUBLE, th,
                                         MB_TAG_CREAT | MB_TAG_DENSE, def_val);
          CHKERR postProc.tag_set_data(th, &nodeHandles[0], nodeHandles.size(),
                                       &getCoordsAtGaussPts()(0, 0));
        }

        auto to_str = [](auto i) {
          return boost::lexical_cast<std::string>(i);
        };
        std::string tag_name_base =
            "Type" + to_str(type) + "Side" + to_str(side);
        cerr << "Tag " << tag_name_base << endl;
        cerr << "Order " << data.getOrder() << endl;

        MatrixDouble trans_base = trans(data.getN());
        if (trans_base.size2() != nodeHandles.size())
          SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "wrong size %d != %d", trans_base.size2(),
                   nodeHandles.size());
        for (size_t rr = 0; rr != trans_base.size1(); ++rr) {
          auto tag_name = tag_name_base + "Base" + to_str(rr);
          Tag th;
          CHKERR postProc.tag_get_handle(tag_name.c_str(), 1, MB_TYPE_DOUBLE,
                                         th, MB_TAG_CREAT | MB_TAG_DENSE,
                                         def_val);
          CHKERR postProc.tag_set_data(th, &nodeHandles[0], nodeHandles.size(),
                                       &trans_base(rr, 0));
        }

        MoFEMFunctionReturn(0);
      }

    public:
      moab::Interface &postProc;
      MapCoords &mapCoords;
      std::vector<EntityHandle> nodeHandles;
    };

    struct MyPrisms : public FatPrismElementForcesAndSourcesCore {

      MyPrisms(MoFEM::Interface &m_field, MatrixDouble &tri_coords)
          : FatPrismElementForcesAndSourcesCore(m_field),
            triCoords(tri_coords) {}

      int getRuleTrianglesOnly(int order) { return -1; };
      int getRuleThroughThickness(int order) { return -1; };

      MoFEMErrorCode setGaussPtsTrianglesOnly(int order_triangles_only) {
        MoFEMFunctionBegin;
        gaussPtsTrianglesOnly.resize(3, triCoords.size1(), false);
        gaussPtsTrianglesOnly.clear();
        for (int gg = 0; gg != triCoords.size1(); ++gg)
          for (int dd = 0; dd != 2; ++dd)
            gaussPtsTrianglesOnly(dd, gg) = triCoords(gg, dd);
        MoFEMFunctionReturn(0);
      }

      MoFEMErrorCode setGaussPtsThroughThickness(int order_thickness) {
        MoFEMFunctionBegin;
        gaussPtsThroughThickness.resize(2, number_of_prisms_layers + 1, false);
        gaussPtsThroughThickness.clear();
        for (int gg = 0; gg != number_of_prisms_layers + 1; ++gg)
          gaussPtsThroughThickness(0, gg) = delta * gg;

        MoFEMFunctionReturn(0);
      }

    private:
      MatrixDouble &triCoords;
    };

    MyPrisms fe1(m_field, tri_coords);
    fe1.getOpPtrVector().push_back(new MyOp(moab, map_coords));
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE1", fe1);

    CHKERR moab.write_file("prism_mesh.vtk", "VTK", "", &meshset, 1);
    CHKERR moab.write_file("one_prism_mesh.vtk", "VTK", "", &one_prism_meshset,
                           1);
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;
}
