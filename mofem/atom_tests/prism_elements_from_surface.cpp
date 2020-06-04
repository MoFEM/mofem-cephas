/** \file prism_elements_from_surface.cpp
  \example prism_elements_from_surface.cpp
  \brief Adding prims on the surface and checking conformity between quads,
  triangles and prism

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

static constexpr int precision_exponent = 5;
static constexpr int number_of_prisms_layers = 18;
static constexpr double delta =
    1. / static_cast<double>(number_of_prisms_layers);
static constexpr std::array<double, 3> d3 = {0, 0, 0};
static constexpr std::array<double, 3> d4 = {0, 0, delta};

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
            CoordsAndHandle, member<CoordsAndHandle, int, &CoordsAndHandle::x>,
            member<CoordsAndHandle, int, &CoordsAndHandle::y>,
            member<CoordsAndHandle, int, &CoordsAndHandle::z>>>

        >>
    MapCoords;

struct PrismOp : public FatPrismElementForcesAndSourcesCore::UserDataOperator {

  PrismOp(moab::Interface &post_proc, MapCoords &map_coords);
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

public:
  moab::Interface &postProc;
  MapCoords &mapCoords;
  std::vector<EntityHandle> nodeHandles;
};

struct PrismFE : public FatPrismElementForcesAndSourcesCore {

  PrismFE(MoFEM::Interface &m_field, MatrixDouble &tri_coords);

  int getRuleTrianglesOnly(int order);
  int getRuleThroughThickness(int order);

  MoFEMErrorCode setGaussPtsTrianglesOnly(int order_triangles_only);

  MoFEMErrorCode setGaussPtsThroughThickness(int order_thickness);

private:
  MatrixDouble &triCoords;
  EntityHandle prims;
};

template <typename OP> struct Op : public OP {

  Op(moab::Interface &post_proc, MapCoords &map_coords, EntityHandle prism);
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

public:
  moab::Interface &postProc;
  MapCoords &mapCoords;
  EntityHandle prism;
  std::vector<EntityHandle> nodeHandles;
};

struct TriFE : public FaceElementForcesAndSourcesCore {

  TriFE(MoFEM::Interface &m_field, MatrixDouble &tri_coords,
        EntityHandle prims);

  int getRule(int order_row, int order_col, int order_data);

  MoFEMErrorCode setGaussPts(int order_row, int order_col, int order_data);

private:
  MatrixDouble &triCoords;
  EntityHandle prism;
};

struct QuadFE : public FaceElementForcesAndSourcesCore {

  QuadFE(MoFEM::Interface &m_field, std::array<Range, 3> &edges_blocks,
         EntityHandle prims);

  int getRule(int order_row, int order_col, int order_data);

  MoFEMErrorCode setGaussPts(int order_row, int order_col, int order_data);

private:
  std::array<Range, 3> &edgeBlocks;
  EntityHandle prism;
};

struct EdgeFE : public EdgeElementForcesAndSourcesCore {

  EdgeFE(MoFEM::Interface &m_field, std::array<Range, 3> &edges_blocks,
         EntityHandle prims);

  int getRule(int order_row, int order_col, int order_data);

  MoFEMErrorCode setGaussPts(int order_row, int order_col, int order_data);

private:
  std::array<Range, 3> &edgeBlocks;
  EntityHandle prism;
};

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
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
              "error -my_file (MESH FILE NEEDED)");

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

    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    std::array<Range, 3> edge_block;
    for (auto b : {1, 2, 3})
      CHKERR m_field.getInterface<MeshsetsManager>()->getEntitiesByDimension(
          b, BLOCKSET, 1, edge_block[b - 1]);

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
    CHKERR m_field.set_field_order(one_prism_meshset, MBEDGE, "FIELD1", 5,
                                   VERY_NOISY);
    CHKERR m_field.set_field_order(one_prism_meshset, MBTRI, "FIELD1", 5);
    CHKERR m_field.set_field_order(one_prism_meshset, MBQUAD, "FIELD1", 5,
                                   VERY_NOISY);
    CHKERR m_field.set_field_order(one_prism_meshset, MBPRISM, "FIELD1", 7,
                                   VERY_NOISY);
    CHKERR m_field.build_fields(VERY_NOISY);

    // FE
    CHKERR m_field.add_finite_element("PRISM");
    CHKERR m_field.add_finite_element("EDGE");
    CHKERR m_field.add_finite_element("TRI");
    CHKERR m_field.add_finite_element("QUAD");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("PRISM", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("PRISM", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("PRISM", "FIELD1");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("EDGE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("EDGE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("EDGE", "FIELD1");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TRI", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("TRI", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TRI", "FIELD1");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("QUAD", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("QUAD", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("QUAD", "FIELD1");

    CHKERR m_field.add_ents_to_finite_element_by_type(one_prism_meshset, MBEDGE,
                                                      "EDGE");
    CHKERR m_field.add_ents_to_finite_element_by_type(one_prism_meshset, MBTRI,
                                                      "TRI");
    CHKERR m_field.add_ents_to_finite_element_by_type(one_prism_meshset, MBQUAD,
                                                      "QUAD");
    CHKERR m_field.add_ents_to_finite_element_by_type(one_prism_meshset,
                                                      MBPRISM, "PRISM");

    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // //build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "PRISM");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "EDGE");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "TRI");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "QUAD");
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

    PrismFE fe_prism(m_field, tri_coords);
    fe_prism.getOpPtrVector().push_back(new PrismOp(moab, map_coords));
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "PRISM", fe_prism);

    EdgeFE fe_edge(m_field, edge_block, one_prism);
    fe_edge.getOpPtrVector().push_back(
        new Op<EdgeElementForcesAndSourcesCoreBase::UserDataOperator>(
            moab, map_coords, one_prism));
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "EDGE", fe_edge);

    TriFE fe_tri(m_field, tri_coords, one_prism);
    fe_tri.getOpPtrVector().push_back(
        new Op<FaceElementForcesAndSourcesCore::UserDataOperator>(
            moab, map_coords, one_prism));
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TRI", fe_tri);

    QuadFE fe_quad(m_field, edge_block, one_prism);
    fe_quad.getOpPtrVector().push_back(
        new Op<FaceElementForcesAndSourcesCoreBase::UserDataOperator>(
            moab, map_coords, one_prism));
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "QUAD", fe_quad);

    CHKERR moab.write_file("prism_mesh.vtk", "VTK", "", &meshset, 1);
    CHKERR moab.write_file("one_prism_mesh.vtk", "VTK", "", &one_prism_meshset,
                           1);
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;
}

PrismOp::PrismOp(moab::Interface &post_proc, MapCoords &map_coords)
    : FatPrismElementForcesAndSourcesCore::UserDataOperator(
          "FIELD1", "FIELD1", ForcesAndSourcesCore::UserDataOperator::OPROW),
      postProc(post_proc), mapCoords(map_coords) {}

MoFEMErrorCode PrismOp::doWork(int side, EntityType type,
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

  const int nb_dofs = data.getIndices().size();
  for (int dd = 0; dd != nb_dofs; ++dd)
    data.getFieldDofs()[dd]->getFieldData() = data.getIndices()[dd];

  if (type == MBVERTEX) {
    const size_t nb_gauss_pts = getGaussPts().size2();
    auto &coords_at_pts = getGaussPts();
    nodeHandles.reserve(nb_gauss_pts);
    nodeHandles.clear();
    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
      auto t = boost::make_tuple(CoordsAndHandle::getArg(coords_at_pts(0, gg)),
                                 CoordsAndHandle::getArg(coords_at_pts(1, gg)),
                                 CoordsAndHandle::getArg(coords_at_pts(2, gg)));

      auto it = mapCoords.find(t);
      if (it != mapCoords.end())
        nodeHandles.emplace_back(it->node);
      else
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Vertex not found");
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
        if (std::abs(node_coords(gg, d) - getCoordsAtGaussPts()(gg, d)) > eps)
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

  auto to_str = [](auto i) { return boost::lexical_cast<std::string>(i); };
  std::string tag_name_base =
      "PrismType" + to_str(type) + "Side" + to_str(side);
  std::cout << tag_name_base << endl;

  MatrixDouble trans_base = trans(data.getN());
  if (trans_base.size2() != nodeHandles.size())
    SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "wrong size %d != %d",
             trans_base.size2(), nodeHandles.size());
  for (size_t rr = 0; rr != trans_base.size1(); ++rr) {
    auto tag_name = tag_name_base + "Base" + to_str(rr);
    Tag th;
    CHKERR postProc.tag_get_handle(tag_name.c_str(), 1, MB_TYPE_DOUBLE, th,
                                   MB_TAG_CREAT | MB_TAG_DENSE, def_val);
    CHKERR postProc.tag_set_data(th, &nodeHandles[0], nodeHandles.size(),
                                 &trans_base(rr, 0));
  }

  MoFEMFunctionReturn(0);
}

PrismFE::PrismFE(MoFEM::Interface &m_field, MatrixDouble &tri_coords)
    : FatPrismElementForcesAndSourcesCore(m_field), triCoords(tri_coords) {}

int PrismFE::getRuleTrianglesOnly(int order) { return -1; };
int PrismFE::getRuleThroughThickness(int order) { return -1; };

MoFEMErrorCode PrismFE::setGaussPtsTrianglesOnly(int order_triangles_only) {
  MoFEMFunctionBegin;
  gaussPtsTrianglesOnly.resize(3, triCoords.size1(), false);
  gaussPtsTrianglesOnly.clear();
  for (int gg = 0; gg != triCoords.size1(); ++gg)
    for (int dd = 0; dd != 2; ++dd)
      gaussPtsTrianglesOnly(dd, gg) = triCoords(gg, dd);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PrismFE::setGaussPtsThroughThickness(int order_thickness) {
  MoFEMFunctionBegin;
  gaussPtsThroughThickness.resize(2, number_of_prisms_layers + 1, false);
  gaussPtsThroughThickness.clear();
  for (int gg = 0; gg != number_of_prisms_layers + 1; ++gg)
    gaussPtsThroughThickness(0, gg) = delta * gg;

  MoFEMFunctionReturn(0);
}

TriFE::TriFE(MoFEM::Interface &m_field, MatrixDouble &tri_coords,
             EntityHandle prism)
    : FaceElementForcesAndSourcesCore(m_field), triCoords(tri_coords),
      prism(prism) {}

int TriFE::getRule(int order_row, int order_col, int order_data) { return -1; }

MoFEMErrorCode TriFE::setGaussPts(int order_row, int order_col,
                                  int order_data) {
  MoFEMFunctionBegin;

  const EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  int side, sense, offset;
  CHKERR mField.get_moab().side_number(prism, ent, side, sense, offset);
  std::array<int, 2> swap = {0, 1};
  if (side == 3)
    swap = std::array<int, 2>{1, 0};

  gaussPts.resize(3, triCoords.size1(), false);
  gaussPts.clear();
  for (int gg = 0; gg != triCoords.size1(); ++gg)
    for (int dd = 0; dd != 2; ++dd)
      gaussPts(dd, gg) = triCoords(gg, swap[dd]);

  MoFEMFunctionReturn(0);
}

QuadFE::QuadFE(MoFEM::Interface &m_field, std::array<Range, 3> &edge_blocks,
               EntityHandle prism)
    : FaceElementForcesAndSourcesCore(m_field), edgeBlocks(edge_blocks),
      prism(prism) {}

int QuadFE::getRule(int order_row, int order_col, int order_data) { return -1; }

MoFEMErrorCode QuadFE::setGaussPts(int order_row, int order_col,
                                   int order_data) {
  MoFEMFunctionBegin;

  const EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  int side, sense, offset;
  CHKERR mField.get_moab().side_number(prism, ent, side, sense, offset);

  Range edge_verts;
  CHKERR mField.get_moab().get_connectivity(edgeBlocks[side], edge_verts);
  MatrixDouble edge_coords(edge_verts.size(), 3);
  CHKERR mField.get_moab().get_coords(edge_verts, &edge_coords(0, 0));

  constexpr double normal[3][2] = {{-1, 0}, {-0.5, 0.5}, {0, 1}};
  constexpr double origin[3][2] = {{1, 0}, {1, 0}, {0, 0}};
  constexpr int swap[3][2] = {{0, 1}, {0, 1}, {1, 0}};
  gaussPts.resize(3, edge_verts.size() * (number_of_prisms_layers + 1), false);
  int gg = 0;
  for (size_t rr = 0; rr != edge_verts.size(); ++rr) {
    const double x = edge_coords(rr, 0) - origin[side][0];
    const double y = edge_coords(rr, 1) - origin[side][1];
    const double edge_dist = x * normal[side][0] + y * normal[side][1];
    for (size_t cc = 0; cc != number_of_prisms_layers + 1; ++cc, ++gg) {
      gaussPts(swap[side][0], gg) = edge_dist;
      gaussPts(swap[side][1], gg) = delta * cc;
    }
  }

  MoFEMFunctionReturn(0);
}

EdgeFE::EdgeFE(MoFEM::Interface &m_field, std::array<Range, 3> &edge_blocks,
               EntityHandle prism)
    : EdgeElementForcesAndSourcesCore(m_field), edgeBlocks(edge_blocks),
      prism(prism) {}

int EdgeFE::getRule(int order_row, int order_col, int order_data) { return -1; }

MoFEMErrorCode EdgeFE::setGaussPts(int order_row, int order_col,
                                   int order_data) {
  MoFEMFunctionBegin;

  const EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  int side, sense, offset;
  CHKERR mField.get_moab().side_number(prism, ent, side, sense, offset);

  if (side >= 3 && side <= 5) {

    gaussPts.resize(2, number_of_prisms_layers + 1, false);
    for (size_t gg = 0; gg != number_of_prisms_layers + 1; ++gg)
      gaussPts(0, gg) = delta * gg;

    // cerr << gaussPts << endl;

  } else {

    constexpr double normal[3][2] = {{1, 0}, {-0.5, 0.5}, {0, 1}};
    constexpr double origin[3][2] = {{0, 1}, {1, 0}, {0, 0}};
    constexpr int side_map[9] = {0, 1, 2, -1, -1, -1, 0, 1, 2};

    int num_nodes;
    const EntityHandle *conn;
    CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
    MatrixDouble coords(num_nodes, 3);
    CHKERR mField.get_moab().get_coords(conn, num_nodes,
                                        &*coords.data().begin());
    side = side_map[side];

    Range edge_verts;
    CHKERR mField.get_moab().get_connectivity(edgeBlocks[side], edge_verts);
    MatrixDouble edge_coords(edge_verts.size(), 3);
    CHKERR mField.get_moab().get_coords(edge_verts, &edge_coords(0, 0));
    gaussPts.resize(2, edge_verts.size(), false);

    for (size_t gg = 0; gg != edge_verts.size(); ++gg) {
      const double x = edge_coords(gg, 0) - origin[side][0];
      const double y = edge_coords(gg, 1) - origin[side][1];
      const double edge_dist = x * normal[side][0] + y * normal[side][1];
      gaussPts(0, gg) = edge_dist;
    }
  }

  MoFEMFunctionReturn(0);
}

template <typename OP>
Op<OP>::Op(moab::Interface &post_proc, MapCoords &map_coords,
           EntityHandle prism)
    : OP("FIELD1", "FIELD1", ForcesAndSourcesCore::UserDataOperator::OPROW),
      postProc(post_proc), mapCoords(map_coords), prism(prism) {}

template <typename OP>
MoFEMErrorCode Op<OP>::doWork(int side, EntityType type,
                              DataForcesAndSourcesCore::EntData &data) {
  constexpr double def_val[] = {0, 0, 0};
  MoFEMFunctionBegin;
  switch (type) {
  case MBVERTEX:
  case MBEDGE:
  case MBTRI:
  case MBQUAD:
    break;
  default:
    MoFEMFunctionReturnHot(0);
  }

  const int nb_dofs = data.getIndices().size();
  for (int dd = 0; dd != nb_dofs; ++dd)
    if (data.getFieldData()[dd] != data.getIndices()[dd]) {
      std::cerr << "Indices: " << data.getIndices() << std::endl;
      std::cerr << "Local indices: " << data.getLocalIndices() << std::endl;
      std::cerr << "Data: " << data.getFieldData() << std::endl;
      SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
               "Indicices/data inconsistency %3.1f != %d",
               data.getFieldData()[dd], data.getIndices()[dd]);
    }

  const EntityHandle fe_ent = OP::getFEEntityHandle();
  const EntityHandle ent = OP::getSideEntity(side, type);
  int side_prism, sense, offset;
  if (type == MBVERTEX) {
    CHKERR postProc.side_number(prism, fe_ent, side_prism, sense, offset);
  } else
    CHKERR postProc.side_number(prism, ent, side_prism, sense, offset);

  if (type == MBVERTEX) {
    auto &coords_at_pts = OP::getCoordsAtGaussPts();
    const size_t nb_gauss_pts = coords_at_pts.size1();

    nodeHandles.reserve(nb_gauss_pts);
    nodeHandles.clear();
    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
      auto t = boost::make_tuple(CoordsAndHandle::getArg(coords_at_pts(gg, 0)),
                                 CoordsAndHandle::getArg(coords_at_pts(gg, 1)),
                                 CoordsAndHandle::getArg(coords_at_pts(gg, 2)));

      auto it = mapCoords.find(t);
      if (it != mapCoords.end())
        nodeHandles.emplace_back(it->node);
      else
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Vertex not found");
    }
  }

  auto to_str = [](auto i) { return boost::lexical_cast<std::string>(i); };
  std::string tag_name_base =
      "FEType" + to_str(OP::getNumeredEntFiniteElementPtr()->getEntType()) +
      "Type" + to_str(type) + "Side" + to_str(side_prism);

  std::string tag_prism_name_base =
      "PrismType" + to_str(type) + "Side" + to_str(side_prism);

  MatrixDouble trans_base = trans(data.getN());
  MatrixDouble prism_base(trans_base.size1(), trans_base.size2());
  if (trans_base.size2() != nodeHandles.size())
    SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "wrong size %d != %d",
             trans_base.size2(), nodeHandles.size());

  for (size_t rr = 0; rr != trans_base.size1(); ++rr) {

    std::string tag_name = tag_name_base + "Base";
    if (type == MBVERTEX) {
      EntityHandle node = data.getFieldDofs()[rr]->getEnt();
      int prism_mode_side;
      CHKERR postProc.side_number(prism, node, prism_mode_side, sense, offset);
      tag_name += to_str(prism_mode_side);
    } else {
      tag_name += to_str(rr);
    }

    std::cout << tag_name << endl;

    Tag th;
    CHKERR postProc.tag_get_handle(tag_name.c_str(), 1, MB_TYPE_DOUBLE, th,
                                   MB_TAG_CREAT | MB_TAG_DENSE, def_val);
    CHKERR postProc.tag_set_data(th, &nodeHandles[0], nodeHandles.size(),
                                 &trans_base(rr, 0));

    if (type != MBVERTEX) {
      auto tag_prism_name = tag_prism_name_base + "Base" + to_str(rr);
      Tag th_prism;
      CHKERR postProc.tag_get_handle(tag_prism_name.c_str(), th_prism);
      CHKERR postProc.tag_get_data(th_prism, &nodeHandles[0],
                                   nodeHandles.size(), &prism_base(rr, 0));
    }
  }

  auto sum_matrix = [](MatrixDouble &m) {
    double s = 0;
    for (unsigned int ii = 0; ii < m.size1(); ii++) {
      for (unsigned int jj = 0; jj < m.size2(); jj++) {
        s += std::abs(m(ii, jj));
      }
    }
    return s;
  };

  if (type != MBVERTEX) {
    prism_base -= trans_base;
    double sum = sum_matrix(prism_base);
    constexpr double eps = 1e-6;

    if (std::abs(sum) > eps)
      cout << "Inconsistent base " << tag_prism_name_base << " "
           << tag_name_base << " sum  " << sum << endl;

    if (std::abs(sum) > eps)
      SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
               "Inconsistent base %s sum %6.4e", tag_prism_name_base.c_str(),
               sum);
  }

  MoFEMFunctionReturn(0);
}
