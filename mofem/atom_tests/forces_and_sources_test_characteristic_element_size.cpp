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

    auto simple = m_field.getInterface<Simple>();

    CHKERR simple->getOptions();
    CHKERR simple->loadFile();

    int order = 1;
    MOFEM_LOG("WORLD", Sev::inform) << "Set order: " << order;
    // CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-DIM", &DIM, PETSC_NULL);

    int base_number = 0;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-base_num", &base_number,
                              PETSC_NULL);

    FieldApproximationBase base = NOBASE;
    if (base_number == 0)
      base = AINSWORTH_LEGENDRE_BASE;
    else if (base_number == 1)
      base = DEMKOWICZ_JACOBI_BASE;

    CHKERR simple->addDomainField("FIELD1", H1, base, SPACE_DIM);
    CHKERR simple->setFieldOrder("FIELD1", order);
    // simple->getAddDomainFE() = true;
    CHKERR simple->setUp();

    struct TestOperator : public DomainEleOp {

      TestOperator() : DomainEleOp("FIELD1", OPROW) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
        MoFEMFunctionBegin;

        constexpr double tol = std::numeric_limits<double>::epsilon();
        const auto type = getNumeredEntFiniteElementPtr()->getEntType();
        const double characteristic_length = getElementCharacteristicLength();
        const double element_measure = getMeasure();
        if (type == MBTRI || type == MBQUAD) {
          if (abs(characteristic_length - sqrt(2.)) > tol)
            SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                     "The characteristic length of face element is wrong, it "
                     "is %e while should be %e",
                     characteristic_length, sqrt(2.));
          if (type == MBTRI)
            if (abs(element_measure - 0.5) > tol)
              SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "The measure of face test element is wrong, it "
                       "is %e while should be %e",
                       characteristic_length, 0.5);

          if (type == MBQUAD)
            if (abs(element_measure - 1.) > tol)
              SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "The measure of face test element is wrong, it "
                       "is %e while should be %e",
                       characteristic_length, 1.);

        } else if (type == MBHEX || type == MBTET) {
          if (abs(characteristic_length - sqrt(3.)) > tol)
            SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                     "The characteristic length of Volume element is wrong, it "
                     "is %e while should be %e",
                     characteristic_length, sqrt(2.));

          if (type == MBHEX)
            if (abs(element_measure - 1) > tol) {
              SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "The measure of MBHEX test element is wrong, it "
                       "is %e while should be %s",
                       characteristic_length, "1");
            }

          if (type == MBTET)
            if (abs(element_measure - 1. / 6.) > tol)
              SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "The measure of MBTET test element is wrong, it "
                       "is %e while should be %e",
                       characteristic_length, 0.5);
        }

        MoFEMFunctionReturn(0);
      }
    };

    auto pipeline_mng = m_field.getInterface<PipelineManager>();

    auto op_bdy_fe = new TestOperator();

    pipeline_mng->getOpDomainRhsPipeline().push_back(op_bdy_fe);
    CHKERR pipeline_mng->loopFiniteElements();
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

  return 0;
}