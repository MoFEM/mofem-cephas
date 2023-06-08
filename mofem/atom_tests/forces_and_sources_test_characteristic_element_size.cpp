/** \file forces_and_sources_test_characteristic_element_size.cpp

  \brief Atom test checking correctness of evaluating characteristic element
  size are calculated on elements

*/

#include <MoFEM.hpp>

using namespace MoFEM;

constexpr int SPACE_DIM =
    EXECUTABLE_DIMENSION; //< Space dimension of problem, mesh

#include <boost/math/quadrature/gauss_kronrod.hpp>
using namespace boost::math::quadrature;

template <int DIM> struct ElementsAndOps {};

template <> struct ElementsAndOps<2> {
  using DomainEle = FaceElementForcesAndSourcesCore;
  using DomainEleOp = DomainEle::UserDataOperator;
};

template <> struct ElementsAndOps<3> {
  using DomainEle = VolumeElementForcesAndSourcesCore;
  using DomainEleOp = DomainEle::UserDataOperator;
};

// using EntData = EntitiesFieldData::EntData;
using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainEleOp = ElementsAndOps<SPACE_DIM>::DomainEleOp;

namespace bio = boost::iostreams;
using bio::stream;
using bio::tee_device;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    //! [Register MoFEM discrete manager in PETSc]
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);
    //! [Register MoFEM discrete manager in PETSc

    //! [Create MoAB]
    moab::Core mb_instance;              ///< mesh database
    moab::Interface &moab = mb_instance; ///< mesh database interface
    //! [Create MoAB]

    //! [Create MoFEM]
    MoFEM::Core core(moab);           ///< finite element database
    MoFEM::Interface &m_field = core; ///< finite element database interface
    //! [Create MoFEM]

    auto generate_adjacencies = [&](EntityHandle &geom_element) {
      MoFEMFunctionBeginHot;

      Range one_geom_ele_range;
      one_geom_ele_range.insert(geom_element);
      Range one_geom_ele_adj_ents_edges;
      CHKERR moab.get_adjacencies(&geom_element, 1, 1, true,
                                  one_geom_ele_adj_ents_edges/*,
                                  moab::Interface::UNION*/);
      if (SPACE_DIM == 3) {
        Range one_geo_ele_adj_ents_faces;
        CHKERR moab.get_adjacencies(&geom_element, 1,  2, true,
                                  one_geo_ele_adj_ents_faces/*,
                                  moab::Interface::UNION*/);
      }

      MoFEMFunctionReturnHot(0);
    };

    int test_num = -1;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-test_num", &test_num,
                              PETSC_NULL);

    switch (test_num) {
    case 0: {
      double tet_coords[] = {-0.5, -0.5, -0.5, 0.5,  -0.5, -0.5,
                             -0.5, 0.5,  -0.5, -0.5, -0.5, 0.5};
      EntityHandle nodes[4];
      for (int nn = 0; nn < 4; nn++) {
        cerr << "tet_coords[3 * nn]" << tet_coords[3 * nn] << "\n";
        CHKERR moab.create_vertex(&tet_coords[3 * nn], nodes[nn]);
      }
      EntityHandle tet;
      CHKERR moab.create_element(MBTET, nodes, 4, tet);
      CHKERR generate_adjacencies(tet);
      break;
    }
    case 1: {
      double hex_coords[] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                             0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
      EntityHandle nodes[8];
      for (int nn = 0; nn < 8; nn++) {
        CHKERR moab.create_vertex(&hex_coords[3 * nn], nodes[nn]);
      }
      EntityHandle hex;
      CHKERR moab.create_element(MBHEX, nodes, 8, hex);
      CHKERR generate_adjacencies(hex);
      break;
    }
    case 2: {
      double tri_coords[] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
      EntityHandle nodes[3];
      for (int nn = 0; nn < 3; nn++) {
        CHKERR moab.create_vertex(&tri_coords[3 * nn], nodes[nn]);
      }
      EntityHandle tri;
      CHKERR moab.create_element(MBTRI, nodes, 3, tri);
      CHKERR generate_adjacencies(tri);
      break;
    }
    case 3: {
      constexpr double base_coords[] = {

          0, 0,
          0,

          1, 0,
          0,

          1, 1,
          0,

          0, 1,
          0

      };
      EntityHandle nodes[4];
      for (int nn = 0; nn < 4; nn++)
        CHKERR moab.create_vertex(&base_coords[3 * nn], nodes[nn]);
      EntityHandle quad;
      CHKERR moab.create_element(MBQUAD, nodes, 4, quad);
      CHKERR generate_adjacencies(quad);
      break;
    }
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
    }

    BitRefLevel bit_level0 = BitRefLevel().set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, SPACE_DIM, bit_level0);

    auto simple = m_field.getInterface<Simple>();
    simple->setDim(SPACE_DIM);

    int order = 1;
    MOFEM_LOG("WORLD", Sev::inform) << "Set order: " << order;

    int base_number = 0;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-base_num", &base_number,
                              PETSC_NULL);

    FieldApproximationBase base = NOBASE;

    if (base_number == 0) {
      base = AINSWORTH_LEGENDRE_BASE;
    } else if (base_number == 1) {
      base = DEMKOWICZ_JACOBI_BASE;
    }

    CHKERR simple->addDomainField("FIELD1", H1, base, SPACE_DIM);
    CHKERR simple->setFieldOrder("FIELD1", order);
    CHKERR simple->setUp();

    struct TestOperator : public DomainEleOp {

      TestOperator(boost::shared_ptr<double> char_length)
          : DomainEleOp("FIELD1", OPROW), charLength(char_length) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
        MoFEMFunctionBegin;

        constexpr double tol = std::numeric_limits<double>::epsilon();
        const auto type = getNumeredEntFiniteElementPtr()->getEntType();
        const double characteristic_length = *charLength;
        const double element_measure = getMeasure();
        if (type == MBTRI || type == MBQUAD) {
          if (abs(characteristic_length - sqrt(2.)) > tol)
            SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                     "The characteristic length of face element is wrong, it "
                     "is %e while should be %e",
                     characteristic_length, sqrt(2.));
          switch (type) {
          case MBTRI:
            if (abs(element_measure - 0.5) > tol)
              SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "The measure of face test element is wrong, it "
                       "is %e while should be %e",
                       characteristic_length, 0.5);
            break;
          case MBQUAD:
            if (abs(element_measure - 1.) > tol)
              SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "The measure of face test element is wrong, it "
                       "is %e while should be %e",
                       characteristic_length, 1.);
            break;
          }
        } else if (type == MBHEX || type == MBTET) {
          if (abs(characteristic_length - sqrt(3.)) > tol)
            SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                     "The characteristic length of Volume element is wrong, it "
                     "is %e while should be %e",
                     characteristic_length, sqrt(2.));
          switch (type) {
          case MBHEX:
            if (abs(element_measure - 1) > tol) {
              SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "The measure of MBHEX test element is wrong, it "
                       "is %e while should be %s",
                       characteristic_length, "1");
            }
            break;
          case MBTET:
            if (abs(element_measure - 1. / 6.) > tol)
              SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "The measure of MBTET test element is wrong, it "
                       "is %e while should be %e",
                       characteristic_length, 0.5);
            break;
          }
        }

        MoFEMFunctionReturn(0);
      }

    protected:
      boost::shared_ptr<double> charLength;
    };

    auto pipeline_mng = m_field.getInterface<PipelineManager>();

    auto char_length = boost::make_shared<double>(0);
    auto op_char_length_eval = new OpCalculateCharLength(m_field, char_length);
    auto op_bdy_fe = new TestOperator(char_length);

    pipeline_mng->getOpDomainRhsPipeline().push_back(op_char_length_eval);
    pipeline_mng->getOpDomainRhsPipeline().push_back(op_bdy_fe);
    CHKERR pipeline_mng->loopFiniteElements();
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

  return 0;
}