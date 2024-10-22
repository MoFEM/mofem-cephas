/**
 * \file plastic.cpp
 * \example plastic.cpp
 *
 * Plasticity in 2d and 3d
 *
 */

/* The above code is a preprocessor directive in C++ that checks if the macro
"EXECUTABLE_DIMENSION" has been defined. If it has not been defined, it replaces
the " */
#ifndef EXECUTABLE_DIMENSION
  #define EXECUTABLE_DIMENSION 3
#endif

// #undef ADD_CONTACT

#include <MoFEM.hpp>
#include <MatrixFunction.hpp>
#include <IntegrationRules.hpp>

using namespace MoFEM;

template <int DIM> struct ElementsAndOps;

template <> struct ElementsAndOps<2> {
  using DomainEle = PipelineManager::FaceEle;
  using BoundaryEle = PipelineManager::EdgeEle;
  using SideEle = PipelineManager::ElementsAndOpsByDim<2>::FaceSideEle;
  static constexpr FieldSpace CONTACT_SPACE = HCURL;
};

template <> struct ElementsAndOps<3> {
  using DomainEle = PipelineManager::VolEle;
  using BoundaryEle = PipelineManager::FaceEle;
  using SideEle = PipelineManager::ElementsAndOpsByDim<3>::FaceSideEle;
  static constexpr FieldSpace CONTACT_SPACE = HDIV;
};

constexpr int SPACE_DIM =
    EXECUTABLE_DIMENSION; //< Space dimension of problem, mesh
constexpr auto size_symm = (SPACE_DIM * (SPACE_DIM + 1)) / 2;

constexpr AssemblyType AT =
    (SCHUR_ASSEMBLE) ? AssemblyType::BLOCK_SCHUR
                     : AssemblyType::PETSC; //< selected assembly type
constexpr IntegrationType IT =
    IntegrationType::GAUSS; //< selected integration type

constexpr FieldSpace ElementsAndOps<2>::CONTACT_SPACE;
constexpr FieldSpace ElementsAndOps<3>::CONTACT_SPACE;
constexpr FieldSpace CONTACT_SPACE = ElementsAndOps<SPACE_DIM>::CONTACT_SPACE;

using EntData = EntitiesFieldData::EntData;
using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainEleOp = DomainEle::UserDataOperator;
using BoundaryEle = ElementsAndOps<SPACE_DIM>::BoundaryEle;
using BoundaryEleOp = BoundaryEle::UserDataOperator;
using PostProcEle = PostProcBrokenMeshInMoab<DomainEle>;
using SkinPostProcEle = PostProcBrokenMeshInMoab<BoundaryEle>;
using SideEle = ElementsAndOps<SPACE_DIM>::SideEle;

inline double iso_hardening_exp(double tau, double b_iso) {
  return std::exp(
      std::max(static_cast<double>(std::numeric_limits<float>::min_exponent10),
               -b_iso * tau));
}

/**
 * Isotropic hardening
 */
inline double iso_hardening(double tau, double H, double Qinf, double b_iso,
                            double sigmaY) {
  return H * tau + Qinf * (1. - iso_hardening_exp(tau, b_iso)) + sigmaY;
}

inline double iso_hardening_dtau(double tau, double H, double Qinf,
                                 double b_iso) {
  auto r = [&](auto tau) {
    return H + Qinf * b_iso * iso_hardening_exp(tau, b_iso);
  };
  constexpr double eps = 1e-12;
  return std::max(r(tau), eps * r(0));
}

/**
 * Kinematic hardening
 */
template <typename T, int DIM>
inline auto
kinematic_hardening(FTensor::Tensor2_symmetric<T, DIM> &t_plastic_strain,
                    double C1_k) {
  FTensor::Index<'i', DIM> i;
  FTensor::Index<'j', DIM> j;
  FTensor::Tensor2_symmetric<double, DIM> t_alpha;
  if (C1_k < std::numeric_limits<double>::epsilon()) {
    t_alpha(i, j) = 0;
    return t_alpha;
  }
  t_alpha(i, j) = C1_k * t_plastic_strain(i, j);
  return t_alpha;
}

template <int DIM>
inline auto kinematic_hardening_dplastic_strain(double C1_k) {
  FTensor::Index<'i', DIM> i;
  FTensor::Index<'j', DIM> j;
  FTensor::Index<'k', DIM> k;
  FTensor::Index<'l', DIM> l;
  FTensor::Ddg<double, DIM, DIM> t_diff;
  constexpr auto t_kd = FTensor::Kronecker_Delta_symmetric<int>();
  t_diff(i, j, k, l) = C1_k * (t_kd(i, k) ^ t_kd(j, l)) / 4.;
  return t_diff;
}

PetscBool is_large_strains = PETSC_TRUE; ///< Large strains
PetscBool set_timer = PETSC_FALSE;       ///< Set timer

double scale = 1.;

double young_modulus = 206913; ///< Young modulus
double poisson_ratio = 0.29;   ///< Poisson ratio
double sigmaY = 450;           ///< Yield stress
double H = 129;                ///< Hardening
double visH = 0;               ///< Viscous hardening
double zeta = 5e-2;            ///< Viscous hardening
double Qinf = 265;             ///< Saturation yield stress
double b_iso = 16.93;          ///< Saturation exponent
double C1_k = 0;               ///< Kinematic hardening

double cn0 = 1;
double cn1 = 1;

int order = 2;             ///< Order displacement
int tau_order = order - 2; ///< Order of tau files
int ep_order = order - 1;  ///< Order of ep files
int geom_order = 2;        ///< Order if fixed.

PetscBool is_quasi_static = PETSC_TRUE;
double rho = 0.0;
double alpha_damping = 0;

#include <HenckyOps.hpp>
#include <PlasticOps.hpp>
#include <PlasticNaturalBCs.hpp>

#ifdef ADD_CONTACT
  #ifdef PYTHON_SDF
    #include <boost/python.hpp>
    #include <boost/python/def.hpp>
    #include <boost/python/numpy.hpp>
namespace bp = boost::python;
namespace np = boost::python::numpy;
  #endif

namespace ContactOps {

double cn_contact = 0.1;

}; // namespace ContactOps

  #include <ContactOps.hpp>
#endif // ADD_CONTACT

using DomainRhsBCs = NaturalBC<DomainEleOp>::Assembly<AT>::LinearForm<IT>;
using OpDomainRhsBCs =
    DomainRhsBCs::OpFlux<PlasticOps::DomainBCs, 1, SPACE_DIM>;
using BoundaryRhsBCs = NaturalBC<BoundaryEleOp>::Assembly<AT>::LinearForm<IT>;
using OpBoundaryRhsBCs =
    BoundaryRhsBCs::OpFlux<PlasticOps::BoundaryBCs, 1, SPACE_DIM>;
using BoundaryLhsBCs = NaturalBC<BoundaryEleOp>::Assembly<AT>::BiLinearForm<IT>;
using OpBoundaryLhsBCs =
    BoundaryLhsBCs::OpFlux<PlasticOps::BoundaryBCs, 1, SPACE_DIM>;

using namespace PlasticOps;
using namespace HenckyOps;
struct Example {

  Example(MoFEM::Interface &m_field) : mField(m_field) {}

  MoFEMErrorCode runProblem();

private:
  MoFEM::Interface &mField;

  MoFEMErrorCode setupProblem();
  MoFEMErrorCode createCommonData();
  MoFEMErrorCode bC();
  MoFEMErrorCode OPs();
  MoFEMErrorCode tsSolve();
  MoFEMErrorCode testOperators();

  boost::shared_ptr<DomainEle> reactionFe;

  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uXScatter;
  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uYScatter;
  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uZScatter;

  struct ScaledTimeScale : public MoFEM::TimeScale {
    using MoFEM::TimeScale::TimeScale;
    double getScale(const double time) {
      return scale * MoFEM::TimeScale::getScale(time);
    };
  };

#ifdef ADD_CONTACT
  #ifdef PYTHON_SDF
  boost::shared_ptr<ContactOps::SDFPython> sdfPythonPtr;
  #endif
#endif // ADD_CONTACT
};

//! [Run problem]
MoFEMErrorCode Example::runProblem() {
  MoFEMFunctionBegin;
  CHKERR createCommonData();
  CHKERR setupProblem();
  CHKERR bC();
  CHKERR OPs();
  PetscBool test_ops = PETSC_FALSE;
  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-test_operators", &test_ops,
                             PETSC_NULL);
  if (test_ops == PETSC_FALSE) {
    CHKERR tsSolve();
  } else {
    CHKERR testOperators();
  }
  MoFEMFunctionReturn(0);
}
//! [Run problem]

//! [Set up problem]
MoFEMErrorCode Example::setupProblem() {
  MoFEMFunctionBegin;
  Simple *simple = mField.getInterface<Simple>();

  Range domain_ents;
  CHKERR mField.get_moab().get_entities_by_dimension(0, SPACE_DIM, domain_ents,
                                                     true);
  auto get_ents_by_dim = [&](const auto dim) {
    if (dim == SPACE_DIM) {
      return domain_ents;
    } else {
      Range ents;
      if (dim == 0)
        CHKERR mField.get_moab().get_connectivity(domain_ents, ents, true);
      else
        CHKERR mField.get_moab().get_entities_by_dimension(0, dim, ents, true);
      return ents;
    }
  };

  auto get_base = [&]() {
    auto domain_ents = get_ents_by_dim(SPACE_DIM);
    if (domain_ents.empty())
      CHK_THROW_MESSAGE(MOFEM_NOT_FOUND, "Empty mesh");
    const auto type = type_from_handle(domain_ents[0]);
    switch (type) {
    case MBQUAD:
      return DEMKOWICZ_JACOBI_BASE;
    case MBHEX:
      return DEMKOWICZ_JACOBI_BASE;
    case MBTRI:
      return AINSWORTH_LEGENDRE_BASE;
    case MBTET:
      return AINSWORTH_LEGENDRE_BASE;
    default:
      CHK_THROW_MESSAGE(MOFEM_NOT_FOUND, "Element type not handled");
    }
    return NOBASE;
  };

  const auto base = get_base();
  MOFEM_LOG("PLASTICITY", Sev::inform)
      << "Base " << ApproximationBaseNames[base];

  CHKERR simple->addDomainField("U", H1, base, SPACE_DIM);
  CHKERR simple->addDomainField("EP", L2, base, size_symm);
  CHKERR simple->addDomainField("TAU", L2, base, 1);
  CHKERR simple->addBoundaryField("U", H1, base, SPACE_DIM);

  CHKERR simple->addDataField("GEOMETRY", H1, base, SPACE_DIM);

  PetscBool order_edge = PETSC_FALSE;
  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-order_edge", &order_edge,
                             PETSC_NULL);
  PetscBool order_face = PETSC_FALSE;
  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-order_face", &order_face,
                             PETSC_NULL);
  PetscBool order_volume = PETSC_FALSE;
  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-order_volume", &order_volume,
                             PETSC_NULL);

  if (order_edge || order_face || order_volume) {

    MOFEM_LOG("PLASTICITY", Sev::inform) << "Order edge " << order_edge
        ? "true"
        : "false";
    MOFEM_LOG("PLASTICITY", Sev::inform) << "Order face " << order_face
        ? "true"
        : "false";
    MOFEM_LOG("PLASTICITY", Sev::inform) << "Order volume " << order_volume
        ? "true"
        : "false";

    auto ents = get_ents_by_dim(0);
    if (order_edge)
      ents.merge(get_ents_by_dim(1));
    if (order_face)
      ents.merge(get_ents_by_dim(2));
    if (order_volume)
      ents.merge(get_ents_by_dim(3));
    CHKERR simple->setFieldOrder("U", order, &ents);
  } else {
    CHKERR simple->setFieldOrder("U", order);
  }
  CHKERR simple->setFieldOrder("EP", ep_order);
  CHKERR simple->setFieldOrder("TAU", tau_order);

  CHKERR simple->setFieldOrder("GEOMETRY", geom_order);

#ifdef ADD_CONTACT
  CHKERR simple->addDomainField("SIGMA", CONTACT_SPACE, DEMKOWICZ_JACOBI_BASE,
                                SPACE_DIM);
  CHKERR simple->addBoundaryField("SIGMA", CONTACT_SPACE, DEMKOWICZ_JACOBI_BASE,
                                  SPACE_DIM);

  auto get_skin = [&]() {
    Range body_ents;
    CHKERR mField.get_moab().get_entities_by_dimension(0, SPACE_DIM, body_ents);
    Skinner skin(&mField.get_moab());
    Range skin_ents;
    CHKERR skin.find_skin(0, body_ents, false, skin_ents);
    return skin_ents;
  };

  auto filter_blocks = [&](auto skin) {
    bool is_contact_block = true;
    Range contact_range;
    for (auto m :
         mField.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

             (boost::format("%s(.*)") % "CONTACT").str()

                 ))

    ) {
      is_contact_block =
          true; ///< blocs interation is collective, so that is set irrespective
                ///< if there are entities in given rank or not in the block
      MOFEM_LOG("CONTACT", Sev::inform)
          << "Find contact block set:  " << m->getName();
      auto meshset = m->getMeshset();
      Range contact_meshset_range;
      CHKERR mField.get_moab().get_entities_by_dimension(
          meshset, SPACE_DIM - 1, contact_meshset_range, true);

      CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(
          contact_meshset_range);
      contact_range.merge(contact_meshset_range);
    }
    if (is_contact_block) {
      MOFEM_LOG("SYNC", Sev::inform)
          << "Nb entities in contact surface: " << contact_range.size();
      MOFEM_LOG_SYNCHRONISE(mField.get_comm());
      skin = intersect(skin, contact_range);
    }
    return skin;
  };

  auto filter_true_skin = [&](auto skin) {
    Range boundary_ents;
    ParallelComm *pcomm =
        ParallelComm::get_pcomm(&mField.get_moab(), MYPCOMM_INDEX);
    CHKERR pcomm->filter_pstatus(skin, PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                 PSTATUS_NOT, -1, &boundary_ents);
    return boundary_ents;
  };

  auto boundary_ents = filter_true_skin(filter_blocks(get_skin()));
  CHKERR simple->setFieldOrder("SIGMA", 0);
  CHKERR simple->setFieldOrder("SIGMA", order - 1, &boundary_ents);
#endif

  CHKERR simple->setUp();
  CHKERR simple->addFieldToEmptyFieldBlocks("U", "TAU");

  auto project_ho_geometry = [&]() {
    Projection10NodeCoordsOnField ent_method(mField, "GEOMETRY");
    return mField.loop_dofs("GEOMETRY", ent_method);
  };
  PetscBool project_geometry = PETSC_TRUE;
  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-project_geometry",
                             &project_geometry, PETSC_NULL);
  if (project_geometry) {
    CHKERR project_ho_geometry();
  }

  MoFEMFunctionReturn(0);
}
//! [Set up problem]

//! [Create common data]
MoFEMErrorCode Example::createCommonData() {
  MoFEMFunctionBegin;

  auto get_command_line_parameters = [&]() {
    MoFEMFunctionBegin;

    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-scale", &scale, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-young_modulus",
                                 &young_modulus, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-poisson_ratio",
                                 &poisson_ratio, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-hardening", &H, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-hardening_viscous", &visH,
                                 PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-yield_stress", &sigmaY,
                                 PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-cn0", &cn0, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-cn1", &cn1, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-zeta", &zeta, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-Qinf", &Qinf, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-b_iso", &b_iso, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-C1_k", &C1_k, PETSC_NULL);
    CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-large_strains",
                               &is_large_strains, PETSC_NULL);
    CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-set_timer", &set_timer,
                               PETSC_NULL);

    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);
    PetscBool tau_order_is_set; ///< true if tau order is set
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-tau_order", &tau_order,
                              &tau_order_is_set);
    PetscBool ep_order_is_set; ///< true if tau order is set
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-ep_order", &ep_order,
                              &ep_order_is_set);
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-geom_order", &geom_order,
                              PETSC_NULL);

    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-rho", &rho, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-alpha_damping",
                                 &alpha_damping, PETSC_NULL);

    MOFEM_LOG("PLASTICITY", Sev::inform) << "Young modulus " << young_modulus;
    MOFEM_LOG("PLASTICITY", Sev::inform) << "Poisson ratio " << poisson_ratio;
    MOFEM_LOG("PLASTICITY", Sev::inform) << "Yield stress " << sigmaY;
    MOFEM_LOG("PLASTICITY", Sev::inform) << "Hardening " << H;
    MOFEM_LOG("PLASTICITY", Sev::inform) << "Viscous hardening " << visH;
    MOFEM_LOG("PLASTICITY", Sev::inform) << "Saturation yield stress " << Qinf;
    MOFEM_LOG("PLASTICITY", Sev::inform) << "Saturation exponent " << b_iso;
    MOFEM_LOG("PLASTICITY", Sev::inform) << "Kinematic hardening " << C1_k;
    MOFEM_LOG("PLASTICITY", Sev::inform) << "cn0 " << cn0;
    MOFEM_LOG("PLASTICITY", Sev::inform) << "cn1 " << cn1;
    MOFEM_LOG("PLASTICITY", Sev::inform) << "zeta " << zeta;

    if (tau_order_is_set == PETSC_FALSE)
      tau_order = order - 2;
    if (ep_order_is_set == PETSC_FALSE)
      ep_order = order - 1;

    MOFEM_LOG("PLASTICITY", Sev::inform) << "Approximation order " << order;
    MOFEM_LOG("PLASTICITY", Sev::inform)
        << "Ep approximation order " << ep_order;
    MOFEM_LOG("PLASTICITY", Sev::inform)
        << "Tau approximation order " << tau_order;
    MOFEM_LOG("PLASTICITY", Sev::inform)
        << "Geometry approximation order " << geom_order;

    MOFEM_LOG("PLASTICITY", Sev::inform) << "Density " << rho;
    MOFEM_LOG("PLASTICITY", Sev::inform) << "alpha_damping " << alpha_damping;

    PetscBool is_scale = PETSC_TRUE;
    CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-is_scale", &is_scale,
                               PETSC_NULL);
    if (is_scale) {
      scale /= young_modulus;
    }

    MOFEM_LOG("PLASTICITY", Sev::inform) << "Scale " << scale;

#ifdef ADD_CONTACT
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-cn_contact",
                                 &ContactOps::cn_contact, PETSC_NULL);
    MOFEM_LOG("CONTACT", Sev::inform)
        << "cn_contact " << ContactOps::cn_contact;
#endif // ADD_CONTACT

    CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-quasi_static",
                               &is_quasi_static, PETSC_NULL);
    MOFEM_LOG("PLASTICITY", Sev::inform)
        << "Is quasi static: " << (is_quasi_static ? "true" : "false");

    MoFEMFunctionReturn(0);
  };

  CHKERR get_command_line_parameters();

#ifdef ADD_CONTACT
  #ifdef PYTHON_SDF
  sdfPythonPtr = boost::make_shared<ContactOps::SDFPython>();
  CHKERR sdfPythonPtr->sdfInit("sdf.py");
  ContactOps::sdfPythonWeakPtr = sdfPythonPtr;
  #endif
#endif // ADD_CONTACT

  MoFEMFunctionReturn(0);
}
//! [Create common data]

//! [Boundary condition]
MoFEMErrorCode Example::bC() {
  MoFEMFunctionBegin;

  auto simple = mField.getInterface<Simple>();
  auto bc_mng = mField.getInterface<BcManager>();

  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "REMOVE_X",
                                           "U", 0, 0);
  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "REMOVE_Y",
                                           "U", 1, 1);
  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "REMOVE_Z",
                                           "U", 2, 2);
  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(),
                                           "REMOVE_ALL", "U", 0, 3);

#ifdef ADD_CONTACT
  for (auto b : {"FIX_X", "REMOVE_X"})
    CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), b,
                                             "SIGMA", 0, 0, false, true);
  for (auto b : {"FIX_Y", "REMOVE_Y"})
    CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), b,
                                             "SIGMA", 1, 1, false, true);
  for (auto b : {"FIX_Z", "REMOVE_Z"})
    CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), b,
                                             "SIGMA", 2, 2, false, true);
  for (auto b : {"FIX_ALL", "REMOVE_ALL"})
    CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), b,
                                             "SIGMA", 0, 3, false, true);
  CHKERR bc_mng->removeBlockDOFsOnEntities(
      simple->getProblemName(), "NO_CONTACT", "SIGMA", 0, 3, false, true);
#endif

  CHKERR bc_mng->pushMarkDOFsOnEntities<DisplacementCubitBcData>(
      simple->getProblemName(), "U");

  auto &bc_map = bc_mng->getBcMapByBlockName();
  for (auto bc : bc_map)
    MOFEM_LOG("PLASTICITY", Sev::verbose) << "Marker " << bc.first;

  MoFEMFunctionReturn(0);
}
//! [Boundary condition]

//! [Push operators to pipeline]
MoFEMErrorCode Example::OPs() {
  MoFEMFunctionBegin;
  auto pip_mng = mField.getInterface<PipelineManager>();

  auto integration_rule_bc = [](int, int, int ao) { return 2 * ao; };

  auto vol_rule = [](int, int, int ao) { return 2 * ao + geom_order - 1; };

  auto add_boundary_ops_lhs_mechanical = [&](auto &pip) {
    MoFEMFunctionBegin;

    CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(pip, {HDIV},
                                                              "GEOMETRY");
    pip.push_back(new OpSetHOWeightsOnSubDim<SPACE_DIM>());

    // Add Natural BCs to LHS
    CHKERR BoundaryLhsBCs::AddFluxToPipeline<OpBoundaryLhsBCs>::add(
        pip, mField, "U", Sev::inform);

#ifdef ADD_CONTACT
    auto simple = mField.getInterface<Simple>();
    CHKERR
    ContactOps::opFactoryBoundaryLhs<SPACE_DIM, AT, GAUSS, BoundaryEleOp>(
        pip, "SIGMA", "U");
    CHKERR
    ContactOps::opFactoryBoundaryToDomainLhs<SPACE_DIM, AT, IT, DomainEle>(
        mField, pip, simple->getDomainFEName(), "SIGMA", "U", "GEOMETRY",
        vol_rule);
#endif // ADD_CONTACT

    MoFEMFunctionReturn(0);
  };

  auto add_boundary_ops_rhs_mechanical = [&](auto &pip) {
    MoFEMFunctionBegin;

    CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(pip, {HDIV},
                                                              "GEOMETRY");
    pip.push_back(new OpSetHOWeightsOnSubDim<SPACE_DIM>());

    // Add Natural BCs to RHS
    CHKERR BoundaryRhsBCs::AddFluxToPipeline<OpBoundaryRhsBCs>::add(
        pip, mField, "U", {boost::make_shared<ScaledTimeScale>()}, Sev::inform);

#ifdef ADD_CONTACT
    CHKERR ContactOps::opFactoryBoundaryRhs<SPACE_DIM, AT, IT, BoundaryEleOp>(
        pip, "SIGMA", "U");
#endif // ADD_CONTACT

    MoFEMFunctionReturn(0);
  };

  auto add_domain_ops_lhs = [this](auto &pip) {
    MoFEMFunctionBegin;
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1, HDIV},
                                                          "GEOMETRY");

    if (is_quasi_static == PETSC_FALSE) {

      //! [Only used for dynamics]
      using OpMass = FormsIntegrators<DomainEleOp>::Assembly<AT>::BiLinearForm<
          GAUSS>::OpMass<1, SPACE_DIM>;
      //! [Only used for dynamics]

      auto get_inertia_and_mass_damping = [this](const double, const double,
                                                 const double) {
        auto *pip = mField.getInterface<PipelineManager>();
        auto &fe_domain_lhs = pip->getDomainLhsFE();
        return (rho / scale) * fe_domain_lhs->ts_aa +
               (alpha_damping / scale) * fe_domain_lhs->ts_a;
      };
      pip.push_back(new OpMass("U", "U", get_inertia_and_mass_damping));
    }

    CHKERR PlasticOps::opFactoryDomainLhs<SPACE_DIM, AT, IT, DomainEleOp>(
        mField, "MAT_PLASTIC", pip, "U", "EP", "TAU");

    MoFEMFunctionReturn(0);
  };

  auto add_domain_ops_rhs = [this](auto &pip) {
    MoFEMFunctionBegin;

    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1, HDIV},
                                                          "GEOMETRY");

    CHKERR DomainRhsBCs::AddFluxToPipeline<OpDomainRhsBCs>::add(
        pip, mField, "U",
        {boost::make_shared<ScaledTimeScale>("body_force_hist.txt")},
        Sev::inform);

    // only in case of dynamics
    if (is_quasi_static == PETSC_FALSE) {

      //! [Only used for dynamics]
      using OpInertiaForce = FormsIntegrators<DomainEleOp>::Assembly<
          AT>::LinearForm<IT>::OpBaseTimesVector<1, SPACE_DIM, 1>;
      //! [Only used for dynamics]

      auto mat_acceleration = boost::make_shared<MatrixDouble>();
      pip.push_back(new OpCalculateVectorFieldValuesDotDot<SPACE_DIM>(
          "U", mat_acceleration));
      pip.push_back(
          new OpInertiaForce("U", mat_acceleration, [](double, double, double) {
            return rho / scale;
          }));
      if (alpha_damping > 0) {
        auto mat_velocity = boost::make_shared<MatrixDouble>();
        pip.push_back(
            new OpCalculateVectorFieldValuesDot<SPACE_DIM>("U", mat_velocity));
        pip.push_back(
            new OpInertiaForce("U", mat_velocity, [](double, double, double) {
              return alpha_damping / scale;
            }));
      }
    }

    CHKERR PlasticOps::opFactoryDomainRhs<SPACE_DIM, AT, IT, DomainEleOp>(
        mField, "MAT_PLASTIC", pip, "U", "EP", "TAU");

#ifdef ADD_CONTACT
    CHKERR ContactOps::opFactoryDomainRhs<SPACE_DIM, AT, IT, DomainEleOp>(
        pip, "SIGMA", "U");
#endif // ADD_CONTACT

    MoFEMFunctionReturn(0);
  };

  CHKERR add_domain_ops_lhs(pip_mng->getOpDomainLhsPipeline());
  CHKERR add_domain_ops_rhs(pip_mng->getOpDomainRhsPipeline());

  // Boundary
  CHKERR add_boundary_ops_lhs_mechanical(pip_mng->getOpBoundaryLhsPipeline());
  CHKERR add_boundary_ops_rhs_mechanical(pip_mng->getOpBoundaryRhsPipeline());

  CHKERR pip_mng->setDomainRhsIntegrationRule(vol_rule);
  CHKERR pip_mng->setDomainLhsIntegrationRule(vol_rule);

  CHKERR pip_mng->setBoundaryLhsIntegrationRule(integration_rule_bc);
  CHKERR pip_mng->setBoundaryRhsIntegrationRule(integration_rule_bc);

  auto create_reaction_pipeline = [&](auto &pip) {
    MoFEMFunctionBegin;
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1},
                                                          "GEOMETRY");
    CHKERR PlasticOps::opFactoryDomainReactions<SPACE_DIM, AT, IT, DomainEleOp>(
        mField, "MAT_PLASTIC", pip, "U", "EP", "TAU");
    MoFEMFunctionReturn(0);
  };

  reactionFe = boost::make_shared<DomainEle>(mField);
  reactionFe->getRuleHook = vol_rule;
  CHKERR create_reaction_pipeline(reactionFe->getOpPtrVector());
  reactionFe->postProcessHook =
      EssentialPreProcReaction<DisplacementCubitBcData>(mField, reactionFe);

  MoFEMFunctionReturn(0);
}
//! [Push operators to pipeline]

//! [Solve]
struct SetUpSchur {

  /**
   * @brief Create data structure for handling Schur complement
   *
   * @param m_field
   * @param sub_dm  Schur complement sub dm
   * @param field_split_it IS of Schur block
   * @param ao_map AO map from sub dm to main problem
   * @return boost::shared_ptr<SetUpSchur>
   */
  static boost::shared_ptr<SetUpSchur> createSetUpSchur(

      MoFEM::Interface &m_field, SmartPetscObj<DM> sub_dm,
      SmartPetscObj<IS> field_split_it, SmartPetscObj<AO> ao_map

  );
  virtual MoFEMErrorCode setUp(TS solver) = 0;

protected:
  SetUpSchur() = default;
};

MoFEMErrorCode Example::tsSolve() {
  MoFEMFunctionBegin;

  Simple *simple = mField.getInterface<Simple>();
  PipelineManager *pip_mng = mField.getInterface<PipelineManager>();
  ISManager *is_manager = mField.getInterface<ISManager>();

  auto snes_ctx_ptr = getDMSnesCtx(simple->getDM());

  auto set_section_monitor = [&](auto solver) {
    MoFEMFunctionBegin;
    SNES snes;
    CHKERR TSGetSNES(solver, &snes);
    CHKERR SNESMonitorSet(snes,
                          (MoFEMErrorCode(*)(SNES, PetscInt, PetscReal,
                                             void *))MoFEMSNESMonitorFields,
                          (void *)(snes_ctx_ptr.get()), nullptr);
    MoFEMFunctionReturn(0);
  };

  auto create_post_process_elements = [&]() {
    auto pp_fe = boost::make_shared<PostProcEle>(mField);

    auto push_vol_ops = [this](auto &pip) {
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1, HDIV},
                                                            "GEOMETRY");

      auto [common_plastic_ptr, common_hencky_ptr] =
          PlasticOps::createCommonPlasticOps<SPACE_DIM, IT, DomainEleOp>(
              mField, "MAT_PLASTIC", pip, "U", "EP", "TAU", 1., Sev::inform);

      if (common_hencky_ptr) {
        if (common_plastic_ptr->mGradPtr != common_hencky_ptr->matGradPtr)
          CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY, "Wrong pointer for grad");
      }

      return std::make_pair(common_plastic_ptr, common_hencky_ptr);
    };

    auto push_vol_post_proc_ops = [this](auto &pp_fe, auto &&p) {
      MoFEMFunctionBegin;

      auto &pip = pp_fe->getOpPtrVector();

      auto [common_plastic_ptr, common_hencky_ptr] = p;

      using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

      auto x_ptr = boost::make_shared<MatrixDouble>();
      pip.push_back(
          new OpCalculateVectorFieldValues<SPACE_DIM>("GEOMETRY", x_ptr));
      auto u_ptr = boost::make_shared<MatrixDouble>();
      pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));

      if (is_large_strains) {

        pip.push_back(

            new OpPPMap(

                pp_fe->getPostProcMesh(), pp_fe->getMapGaussPts(),

                {{"PLASTIC_SURFACE",
                  common_plastic_ptr->getPlasticSurfacePtr()},
                 {"PLASTIC_MULTIPLIER",
                  common_plastic_ptr->getPlasticTauPtr()}},

                {{"U", u_ptr}, {"GEOMETRY", x_ptr}},

                {{"GRAD", common_hencky_ptr->matGradPtr},
                 {"FIRST_PIOLA", common_hencky_ptr->getMatFirstPiolaStress()}},

                {{"PLASTIC_STRAIN", common_plastic_ptr->getPlasticStrainPtr()},
                 {"PLASTIC_FLOW", common_plastic_ptr->getPlasticFlowPtr()}}

                )

        );

      } else {

        pip.push_back(

            new OpPPMap(

                pp_fe->getPostProcMesh(), pp_fe->getMapGaussPts(),

                {{"PLASTIC_SURFACE",
                  common_plastic_ptr->getPlasticSurfacePtr()},
                 {"PLASTIC_MULTIPLIER",
                  common_plastic_ptr->getPlasticTauPtr()}},

                {{"U", u_ptr}, {"GEOMETRY", x_ptr}},

                {},

                {{"STRAIN", common_plastic_ptr->mStrainPtr},
                 {"STRESS", common_plastic_ptr->mStressPtr},
                 {"PLASTIC_STRAIN", common_plastic_ptr->getPlasticStrainPtr()},
                 {"PLASTIC_FLOW", common_plastic_ptr->getPlasticFlowPtr()}}

                )

        );
      }

      MoFEMFunctionReturn(0);
    };

    auto vol_post_proc = [this, push_vol_post_proc_ops, push_vol_ops]() {
      PetscBool post_proc_vol = PETSC_FALSE;
      CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-post_proc_vol",
                                 &post_proc_vol, PETSC_NULL);
      if (post_proc_vol == PETSC_FALSE)
        return boost::shared_ptr<PostProcEle>();
      auto pp_fe = boost::make_shared<PostProcEle>(mField);
      CHK_MOAB_THROW(
          push_vol_post_proc_ops(pp_fe, push_vol_ops(pp_fe->getOpPtrVector())),
          "push_vol_post_proc_ops");
      return pp_fe;
    };

    auto skin_post_proc = [this, push_vol_post_proc_ops, push_vol_ops]() {
      PetscBool post_proc_skin = PETSC_TRUE;
      CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-post_proc_skin",
                                 &post_proc_skin, PETSC_NULL);
      if (post_proc_skin == PETSC_FALSE)
        return boost::shared_ptr<SkinPostProcEle>();

      auto simple = mField.getInterface<Simple>();
      auto pp_fe = boost::make_shared<SkinPostProcEle>(mField);
      auto op_side = new OpLoopSide<SideEle>(mField, simple->getDomainFEName(),
                                             SPACE_DIM, Sev::verbose);
      pp_fe->getOpPtrVector().push_back(op_side);
      CHK_MOAB_THROW(push_vol_post_proc_ops(
                         pp_fe, push_vol_ops(op_side->getOpPtrVector())),
                     "push_vol_post_proc_ops");
      return pp_fe;
    };

    return std::make_pair(vol_post_proc(), skin_post_proc());
  };

  auto scatter_create = [&](auto D, auto coeff) {
    SmartPetscObj<IS> is;
    CHKERR is_manager->isCreateProblemFieldAndRank(simple->getProblemName(),
                                                   ROW, "U", coeff, coeff, is);
    int loc_size;
    CHKERR ISGetLocalSize(is, &loc_size);
    Vec v;
    CHKERR VecCreateMPI(mField.get_comm(), loc_size, PETSC_DETERMINE, &v);
    VecScatter scatter;
    CHKERR VecScatterCreate(D, is, v, PETSC_NULL, &scatter);
    return std::make_tuple(SmartPetscObj<Vec>(v),
                           SmartPetscObj<VecScatter>(scatter));
  };

  auto set_time_monitor = [&](auto dm, auto solver) {
    MoFEMFunctionBegin;
    boost::shared_ptr<Monitor> monitor_ptr(
        new Monitor(dm, create_post_process_elements(), reactionFe, uXScatter,
                    uYScatter, uZScatter));
    boost::shared_ptr<ForcesAndSourcesCore> null;
    CHKERR DMMoFEMTSSetMonitor(dm, solver, simple->getDomainFEName(),
                               monitor_ptr, null, null);
    MoFEMFunctionReturn(0);
  };

  auto set_schur_pc = [&](auto solver,
                          boost::shared_ptr<SetUpSchur> &schur_ptr) {
    MoFEMFunctionBeginHot;

    auto name_prb = simple->getProblemName();

    // create sub dm for Schur complement
    auto create_schur_dm = [&](SmartPetscObj<DM> base_dm,
                               SmartPetscObj<DM> &dm_sub) {
      MoFEMFunctionBegin;
      dm_sub = createDM(mField.get_comm(), "DMMOFEM");
      CHKERR DMMoFEMCreateSubDM(dm_sub, base_dm, "SCHUR");
      CHKERR DMMoFEMSetSquareProblem(dm_sub, PETSC_TRUE);
      CHKERR DMMoFEMAddElement(dm_sub, simple->getDomainFEName());
      CHKERR DMMoFEMAddElement(dm_sub, simple->getBoundaryFEName());
      for (auto f : {"U"}) {
        CHKERR DMMoFEMAddSubFieldRow(dm_sub, f);
        CHKERR DMMoFEMAddSubFieldCol(dm_sub, f);
      }
      CHKERR DMSetUp(dm_sub);

      MoFEMFunctionReturn(0);
    };

    auto create_block_dm = [&](SmartPetscObj<DM> base_dm,
                               SmartPetscObj<DM> &dm_sub) {
      MoFEMFunctionBegin;
      dm_sub = createDM(mField.get_comm(), "DMMOFEM");
      CHKERR DMMoFEMCreateSubDM(dm_sub, base_dm, "BLOCK");
      CHKERR DMMoFEMSetSquareProblem(dm_sub, PETSC_TRUE);
      CHKERR DMMoFEMAddElement(dm_sub, simple->getDomainFEName());
      CHKERR DMMoFEMAddElement(dm_sub, simple->getBoundaryFEName());
#ifdef ADD_CONTACT
      for (auto f : {"SIGMA", "EP", "TAU"}) {
        CHKERR DMMoFEMAddSubFieldRow(dm_sub, f);
        CHKERR DMMoFEMAddSubFieldCol(dm_sub, f);
      }
#else
      for (auto f : {"EP", "TAU"}) {
        CHKERR DMMoFEMAddSubFieldRow(dm_sub, f);
        CHKERR DMMoFEMAddSubFieldCol(dm_sub, f);
      }
#endif
      CHKERR DMSetUp(dm_sub);
      MoFEMFunctionReturn(0);
    };

    // Create nested (sub BC) Schur DM
    if constexpr (AT == AssemblyType::BLOCK_SCHUR) {

      SmartPetscObj<DM> dm_schur;
      CHKERR create_schur_dm(simple->getDM(), dm_schur);
      SmartPetscObj<DM> dm_block;
      CHKERR create_block_dm(simple->getDM(), dm_block);

#ifdef ADD_CONTACT

      auto get_nested_mat_data = [&](auto schur_dm, auto block_dm) {
        auto block_mat_data = createBlockMatStructure(
            simple->getDM(),

            {

                {simple->getDomainFEName(),

                 {{"U", "U"},
                  {"SIGMA", "SIGMA"},
                  {"U", "SIGMA"},
                  {"SIGMA", "U"},
                  {"EP", "EP"},
                  {"TAU", "TAU"},
                  {"U", "EP"},
                  {"EP", "U"},
                  {"EP", "TAU"},
                  {"TAU", "EP"},
                  {"TAU", "U"}

                 }},

                {simple->getBoundaryFEName(),

                 {{"SIGMA", "SIGMA"}, {"U", "SIGMA"}, {"SIGMA", "U"}

                 }}

            }

        );

        return createSchurNestedMatrixStruture(

            {dm_schur, dm_block}, block_mat_data,

            {"SIGMA", "EP", "TAU"}, {nullptr, nullptr, nullptr}, true

        );
      };

#else

      auto get_nested_mat_data = [&](auto schur_dm, auto block_dm) {
        auto block_mat_data =
            createBlockMatStructure(simple->getDM(),

                                    {{simple->getDomainFEName(),

                                      {{"U", "U"},
                                       {"EP", "EP"},
                                       {"TAU", "TAU"},
                                       {"U", "EP"},
                                       {"EP", "U"},
                                       {"EP", "TAU"},
                                       {"TAU", "U"},
                                       {"TAU", "EP"}

                                      }}}

            );

        return createSchurNestedMatrixStruture(

            {dm_schur, dm_block}, block_mat_data,

            {"EP", "TAU"}, {nullptr, nullptr}, false

        );
      };

#endif

      auto nested_mat_data = get_nested_mat_data(dm_schur, dm_block);
      CHKERR DMMoFEMSetNestSchurData(simple->getDM(), nested_mat_data);

      auto block_is = getDMSubData(dm_block)->getSmartRowIs();
      auto ao_schur = getDMSubData(dm_schur)->getSmartRowMap();

      // Indices has to be map fro very to level, while assembling Schur
      // complement.
      schur_ptr =
          SetUpSchur::createSetUpSchur(mField, dm_schur, block_is, ao_schur);
      CHKERR schur_ptr->setUp(solver);
    }

    MoFEMFunctionReturnHot(0);
  };

  auto dm = simple->getDM();
  auto D = createDMVector(dm);
  auto DD = vectorDuplicate(D);
  uXScatter = scatter_create(D, 0);
  uYScatter = scatter_create(D, 1);
  if constexpr (SPACE_DIM == 3)
    uZScatter = scatter_create(D, 2);

  auto create_solver = [pip_mng]() {
    if (is_quasi_static == PETSC_TRUE)
      return pip_mng->createTSIM();
    else
      return pip_mng->createTSIM2();
  };

  auto solver = create_solver();

  auto active_pre_lhs = []() {
    MoFEMFunctionBegin;
    std::fill(PlasticOps::CommonData::activityData.begin(),
              PlasticOps::CommonData::activityData.end(), 0);
    MoFEMFunctionReturn(0);
  };

  auto active_post_lhs = [&]() {
    MoFEMFunctionBegin;
    auto get_iter = [&]() {
      SNES snes;
      CHK_THROW_MESSAGE(TSGetSNES(solver, &snes), "Can not get SNES");
      int iter;
      CHK_THROW_MESSAGE(SNESGetIterationNumber(snes, &iter),
                        "Can not get iter");
      return iter;
    };

    auto iter = get_iter();
    if (iter >= 0) {

      std::array<int, 5> activity_data;
      std::fill(activity_data.begin(), activity_data.end(), 0);
      MPI_Allreduce(PlasticOps::CommonData::activityData.data(),
                    activity_data.data(), activity_data.size(), MPI_INT,
                    MPI_SUM, mField.get_comm());

      int &active_points = activity_data[0];
      int &avtive_full_elems = activity_data[1];
      int &avtive_elems = activity_data[2];
      int &nb_points = activity_data[3];
      int &nb_elements = activity_data[4];

      if (nb_points) {

        double proc_nb_points =
            100 * static_cast<double>(active_points) / nb_points;
        double proc_nb_active =
            100 * static_cast<double>(avtive_elems) / nb_elements;
        double proc_nb_full_active = 100;
        if (avtive_elems)
          proc_nb_full_active =
              100 * static_cast<double>(avtive_full_elems) / avtive_elems;

        MOFEM_LOG_C("PLASTICITY", Sev::inform,
                    "Iter %d nb pts %d nb active pts %d (%3.3f\%) nb active "
                    "elements %d "
                    "(%3.3f\%) nb full active elems %d (%3.3f\%)",
                    iter, nb_points, active_points, proc_nb_points,
                    avtive_elems, proc_nb_active, avtive_full_elems,
                    proc_nb_full_active, iter);
      }
    }

    MoFEMFunctionReturn(0);
  };

  auto add_active_dofs_elem = [&](auto dm) {
    MoFEMFunctionBegin;
    auto fe_pre_proc = boost::make_shared<FEMethod>();
    fe_pre_proc->preProcessHook = active_pre_lhs;
    auto fe_post_proc = boost::make_shared<FEMethod>();
    fe_post_proc->postProcessHook = active_post_lhs;
    auto ts_ctx_ptr = getDMTsCtx(dm);
    ts_ctx_ptr->getPreProcessIJacobian().push_front(fe_pre_proc);
    ts_ctx_ptr->getPostProcessIJacobian().push_back(fe_post_proc);
    MoFEMFunctionReturn(0);
  };

  auto set_essential_bc = [&](auto dm, auto solver) {
    MoFEMFunctionBegin;
    // This is low level pushing finite elements (pipelines) to solver

    auto pre_proc_ptr = boost::make_shared<FEMethod>();
    auto post_proc_rhs_ptr = boost::make_shared<FEMethod>();
    auto post_proc_lhs_ptr = boost::make_shared<FEMethod>();

    // Add boundary condition scaling
    auto disp_time_scale = boost::make_shared<TimeScale>();

    auto get_bc_hook_rhs = [this, pre_proc_ptr, disp_time_scale]() {
      EssentialPreProc<DisplacementCubitBcData> hook(mField, pre_proc_ptr,
                                                     {disp_time_scale}, false);
      return hook;
    };
    pre_proc_ptr->preProcessHook = get_bc_hook_rhs();

    auto get_post_proc_hook_rhs = [this, post_proc_rhs_ptr]() {
      MoFEMFunctionBegin;
      CHKERR EssentialPreProcReaction<DisplacementCubitBcData>(
          mField, post_proc_rhs_ptr, nullptr, Sev::verbose)();
      CHKERR EssentialPostProcRhs<DisplacementCubitBcData>(
          mField, post_proc_rhs_ptr, 1.)();
      MoFEMFunctionReturn(0);
    };
    auto get_post_proc_hook_lhs = [this, post_proc_lhs_ptr]() {
      return EssentialPostProcLhs<DisplacementCubitBcData>(
          mField, post_proc_lhs_ptr, 1.);
    };
    post_proc_rhs_ptr->postProcessHook = get_post_proc_hook_rhs;

    auto ts_ctx_ptr = getDMTsCtx(dm);
    ts_ctx_ptr->getPreProcessIFunction().push_front(pre_proc_ptr);
    ts_ctx_ptr->getPreProcessIJacobian().push_front(pre_proc_ptr);
    ts_ctx_ptr->getPostProcessIFunction().push_back(post_proc_rhs_ptr);
    post_proc_lhs_ptr->postProcessHook = get_post_proc_hook_lhs();
    ts_ctx_ptr->getPostProcessIJacobian().push_back(post_proc_lhs_ptr);

    MoFEMFunctionReturn(0);
  };

  if (is_quasi_static == PETSC_TRUE) {
    CHKERR TSSetSolution(solver, D);
  } else {
    CHKERR TS2SetSolution(solver, D, DD);
  }

  CHKERR set_section_monitor(solver);
  CHKERR set_time_monitor(dm, solver);
  CHKERR TSSetFromOptions(solver);
  CHKERR TSSetUp(solver);

  CHKERR add_active_dofs_elem(dm);
  boost::shared_ptr<SetUpSchur> schur_ptr;
  CHKERR set_schur_pc(solver, schur_ptr);
  CHKERR set_essential_bc(dm, solver);

  MOFEM_LOG_CHANNEL("TIMER");
  MOFEM_LOG_TAG("TIMER", "timer");
  if (set_timer)
    BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
  MOFEM_LOG("TIMER", Sev::verbose) << "TSSetUp";
  CHKERR TSSetUp(solver);
  MOFEM_LOG("TIMER", Sev::verbose) << "TSSetUp <= done";
  MOFEM_LOG("TIMER", Sev::verbose) << "TSSolve";
  CHKERR TSSolve(solver, NULL);
  MOFEM_LOG("TIMER", Sev::verbose) << "TSSolve <= done";

  MoFEMFunctionReturn(0);
}
//! [Solve]

//! [TestOperators]
MoFEMErrorCode Example::testOperators() {
  MoFEMFunctionBegin;

  // get operators tester
  auto simple = mField.getInterface<Simple>();
  auto opt = mField.getInterface<OperatorsTester>(); // get interface to
                                                     // OperatorsTester
  auto pip = mField.getInterface<PipelineManager>(); // get interface to
                                                     // pipeline manager

  constexpr double eps = 1e-9;

  auto x = opt->setRandomFields(simple->getDM(), {

                                                     {"U", {-1e-6, 1e-6}},

                                                     {"EP", {-1e-6, 1e-6}},

                                                     {"TAU", {0, 1e-4}}

                                                 });

  auto dot_x_plastic_active =
      opt->setRandomFields(simple->getDM(), {

                                                {"U", {-1, 1}},

                                                {"EP", {-1, 1}},

                                                {"TAU", {0.1, 1}}

                                            });
  auto diff_x_plastic_active =
      opt->setRandomFields(simple->getDM(), {

                                                {"U", {-1, 1}},

                                                {"EP", {-1, 1}},

                                                {"TAU", {-1, 1}}

                                            });

  auto dot_x_elastic =
      opt->setRandomFields(simple->getDM(), {

                                                {"U", {-1, 1}},

                                                {"EP", {-1, 1}},

                                                {"TAU", {-1, -0.1}}

                                            });
  auto diff_x_elastic =
      opt->setRandomFields(simple->getDM(), {

                                                {"U", {-1, 1}},

                                                {"EP", {-1, 1}},

                                                {"TAU", {-1, 1}}

                                            });

  auto test_domain_ops = [&](auto fe_name, auto lhs_pipeline, auto rhs_pipeline,
                             auto dot_x, auto diff_x) {
    MoFEMFunctionBegin;

    auto diff_res = opt->checkCentralFiniteDifference(
        simple->getDM(), fe_name, rhs_pipeline, lhs_pipeline, x, dot_x,
        SmartPetscObj<Vec>(), diff_x, 0, 0.5, eps);

    // Calculate norm of difference between directional derivative calculated 
    // from finite difference, and tangent matrix.
    double fnorm;
    CHKERR VecNorm(diff_res, NORM_2, &fnorm);
    MOFEM_LOG_C("PLASTICITY", Sev::inform,
                "Test consistency of tangent matrix %3.4e", fnorm);

    constexpr double err = 1e-5;
    if (fnorm > err)
      SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
               "Norm of directional derivative too large err = %3.4e", fnorm);

    MoFEMFunctionReturn(0);
  };

  MOFEM_LOG("PLASTICITY", Sev::inform) << "Elastic active";
  CHKERR test_domain_ops(simple->getDomainFEName(), pip->getDomainLhsFE(),
                         pip->getDomainRhsFE(), dot_x_elastic, diff_x_elastic);

  MOFEM_LOG("PLASTICITY", Sev::inform) << "Plastic active";
  CHKERR test_domain_ops(simple->getDomainFEName(), pip->getDomainLhsFE(),
                         pip->getDomainRhsFE(), dot_x_plastic_active,
                         diff_x_plastic_active);



  MoFEMFunctionReturn(0);
};

//! [TestOperators]

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

#ifdef ADD_CONTACT
  #ifdef PYTHON_SDF
  Py_Initialize();
  np::initialize();
  #endif
#endif // ADD_CONTACT

  // Initialisation of MoFEM/PETSc and MOAB data structures
  const char param_file[] = "param_file.petsc";
  MoFEM::Core::Initialize(&argc, &argv, param_file, help);

  // Add logging channel for example
  auto core_log = logging::core::get();
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmWorld(), "PLASTICITY"));
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmWorld(), "TIMER"));
  LogManager::setLog("PLASTICITY");
  MOFEM_LOG_TAG("PLASTICITY", "Plasticity");

#ifdef ADD_CONTACT
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmWorld(), "CONTACT"));
  LogManager::setLog("CONTACT");
  MOFEM_LOG_TAG("CONTACT", "Contact");
#endif // ADD_CONTACT

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

    //! [Load mesh]
    Simple *simple = m_field.getInterface<Simple>();
    CHKERR simple->getOptions();
    CHKERR simple->loadFile();
    //! [Load mesh]

    //! [Example]
    Example ex(m_field);
    CHKERR ex.runProblem();
    //! [Example]
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

#ifdef ADD_CONTACT
  #ifdef PYTHON_SDF
  if (Py_FinalizeEx() < 0) {
    exit(120);
  }
  #endif
#endif // ADD_CONTACT

  return 0;
}

struct SetUpSchurImpl : public SetUpSchur {

  SetUpSchurImpl(MoFEM::Interface &m_field, SmartPetscObj<DM> sub_dm,
                 SmartPetscObj<IS> field_split_is, SmartPetscObj<AO> ao_up)
      : SetUpSchur(), mField(m_field), subDM(sub_dm),
        fieldSplitIS(field_split_is), aoSchur(ao_up) {
    if (S) {
      CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                        "Is expected that schur matrix is not "
                        "allocated. This is "
                        "possible only is PC is set up twice");
    }
  }
  virtual ~SetUpSchurImpl() { S.reset(); }

  MoFEMErrorCode setUp(TS solver);
  MoFEMErrorCode preProc();
  MoFEMErrorCode postProc();

private:
  SmartPetscObj<Mat> S;

  MoFEM::Interface &mField;
  SmartPetscObj<DM> subDM;        ///< field split sub dm
  SmartPetscObj<IS> fieldSplitIS; ///< IS for split Schur block
  SmartPetscObj<AO> aoSchur;      ///> main DM to subDM
};

MoFEMErrorCode SetUpSchurImpl::setUp(TS solver) {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  auto pip_mng = mField.getInterface<PipelineManager>();

  SNES snes;
  CHKERR TSGetSNES(solver, &snes);
  KSP ksp;
  CHKERR SNESGetKSP(snes, &ksp);
  CHKERR KSPSetFromOptions(ksp);

  PC pc;
  CHKERR KSPGetPC(ksp, &pc);
  PetscBool is_pcfs = PETSC_FALSE;
  PetscObjectTypeCompare((PetscObject)pc, PCFIELDSPLIT, &is_pcfs);
  if (is_pcfs) {
    if (S) {
      CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                        "Is expected that schur matrix is not "
                        "allocated. This is "
                        "possible only is PC is set up twice");
    }

    S = createDMMatrix(subDM);
    CHKERR MatSetBlockSize(S, SPACE_DIM);

    // Set DM to use shell block matrix
    DM solver_dm;
    CHKERR TSGetDM(solver, &solver_dm);
    CHKERR DMSetMatType(solver_dm, MATSHELL);

    auto ts_ctx_ptr = getDMTsCtx(solver_dm);
    auto A = createDMBlockMat(simple->getDM());
    auto P = createDMNestSchurMat(simple->getDM());

    if (is_quasi_static == PETSC_TRUE) {
      auto swap_assemble = [](TS ts, PetscReal t, Vec u, Vec u_t, PetscReal a,
                              Mat A, Mat B, void *ctx) {
        return TsSetIJacobian(ts, t, u, u_t, a, B, A, ctx);
      };
      CHKERR TSSetIJacobian(solver, A, P, swap_assemble, ts_ctx_ptr.get());
    } else {
      auto swap_assemble = [](TS ts, PetscReal t, Vec u, Vec u_t, Vec utt,
                              PetscReal a, PetscReal aa, Mat A, Mat B,
                              void *ctx) {
        return TsSetI2Jacobian(ts, t, u, u_t, utt, a, aa, B, A, ctx);
      };
      CHKERR TSSetI2Jacobian(solver, A, P, swap_assemble, ts_ctx_ptr.get());
    }
    CHKERR KSPSetOperators(ksp, A, P);

    auto set_ops = [&]() {
      MoFEMFunctionBegin;
      auto pip_mng = mField.getInterface<PipelineManager>();

#ifndef ADD_CONTACT
      // Boundary
      pip_mng->getOpBoundaryLhsPipeline().push_front(
          createOpSchurAssembleBegin());
      pip_mng->getOpBoundaryLhsPipeline().push_back(createOpSchurAssembleEnd(

          {"EP", "TAU"}, {nullptr, nullptr}, aoSchur, S, false, false

          ));
      // Domain
      pip_mng->getOpDomainLhsPipeline().push_front(
          createOpSchurAssembleBegin());
      pip_mng->getOpDomainLhsPipeline().push_back(createOpSchurAssembleEnd(

          {"EP", "TAU"}, {nullptr, nullptr}, aoSchur, S, false, false

          ));
#else

      double eps_stab = 1e-4;
      CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-eps_stab", &eps_stab,
                                   PETSC_NULL);

      using B = FormsIntegrators<BoundaryEleOp>::Assembly<
          BLOCK_PRECONDITIONER_SCHUR>::BiLinearForm<IT>;
      using OpMassStab = B::OpMass<3, SPACE_DIM * SPACE_DIM>;

      // Boundary
      pip_mng->getOpBoundaryLhsPipeline().push_front(
          createOpSchurAssembleBegin());
      pip_mng->getOpBoundaryLhsPipeline().push_back(
          new OpMassStab("SIGMA", "SIGMA", [eps_stab](double, double, double) {
            return eps_stab;
          }));
      pip_mng->getOpBoundaryLhsPipeline().push_back(createOpSchurAssembleEnd(

          {"SIGMA", "EP", "TAU"}, {nullptr, nullptr, nullptr}, aoSchur, S,
          false, false

          ));
      // Domain
      pip_mng->getOpDomainLhsPipeline().push_front(
          createOpSchurAssembleBegin());
      pip_mng->getOpDomainLhsPipeline().push_back(createOpSchurAssembleEnd(

          {"SIGMA", "EP", "TAU"}, {nullptr, nullptr, nullptr}, aoSchur, S,
          false, false

          ));
#endif // ADD_CONTACT
      MoFEMFunctionReturn(0);
    };

    auto set_assemble_elems = [&]() {
      MoFEMFunctionBegin;
      auto schur_asmb_pre_proc = boost::make_shared<FEMethod>();
      schur_asmb_pre_proc->preProcessHook = [this]() {
        MoFEMFunctionBegin;
        CHKERR MatZeroEntries(S);
        MOFEM_LOG("TIMER", Sev::verbose) << "Lhs Assemble Begin";
        MoFEMFunctionReturn(0);
      };
      auto schur_asmb_post_proc = boost::make_shared<FEMethod>();

      schur_asmb_post_proc->postProcessHook = [this, schur_asmb_post_proc]() {
        MoFEMFunctionBegin;
        MOFEM_LOG("TIMER", Sev::verbose) << "Lhs Assemble End";

        // Apply essential constrains to Schur complement
        CHKERR MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
        CHKERR MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
        CHKERR EssentialPostProcLhs<DisplacementCubitBcData>(
            mField, schur_asmb_post_proc, 1, S, aoSchur)();

        MoFEMFunctionReturn(0);
      };
      auto ts_ctx_ptr = getDMTsCtx(simple->getDM());
      ts_ctx_ptr->getPreProcessIJacobian().push_front(schur_asmb_pre_proc);
      ts_ctx_ptr->getPostProcessIJacobian().push_front(schur_asmb_post_proc);
      MoFEMFunctionReturn(0);
    };

    auto set_pc = [&]() {
      MoFEMFunctionBegin;
      CHKERR PCFieldSplitSetIS(pc, NULL, fieldSplitIS);
      CHKERR PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, S);
      MoFEMFunctionReturn(0);
    };

    auto set_diagonal_pc = [&]() {
      MoFEMFunctionBegin;
      KSP *subksp;
      CHKERR PCFieldSplitSchurGetSubKSP(pc, PETSC_NULL, &subksp);
      auto get_pc = [](auto ksp) {
        PC pc_raw;
        CHKERR KSPGetPC(ksp, &pc_raw);
        return SmartPetscObj<PC>(pc_raw,
                                 true); // bump reference
      };
      CHKERR setSchurA00MatSolvePC(get_pc(subksp[0]));
      CHKERR PetscFree(subksp);
      MoFEMFunctionReturn(0);
    };

    CHKERR set_ops();
    CHKERR set_pc();
    CHKERR set_assemble_elems();

    CHKERR TSSetUp(solver);
    CHKERR KSPSetUp(ksp);
    CHKERR set_diagonal_pc();

  } else {
    pip_mng->getOpBoundaryLhsPipeline().push_front(
        createOpSchurAssembleBegin());
    pip_mng->getOpBoundaryLhsPipeline().push_back(
        createOpSchurAssembleEnd({}, {}));
    pip_mng->getOpDomainLhsPipeline().push_front(createOpSchurAssembleBegin());
    pip_mng->getOpDomainLhsPipeline().push_back(
        createOpSchurAssembleEnd({}, {}));
  }

  // fieldSplitIS.reset();
  // aoSchur.reset();
  MoFEMFunctionReturn(0);
}

boost::shared_ptr<SetUpSchur>
SetUpSchur::createSetUpSchur(MoFEM::Interface &m_field,
                             SmartPetscObj<DM> sub_dm, SmartPetscObj<IS> is_sub,
                             SmartPetscObj<AO> ao_up) {
  return boost::shared_ptr<SetUpSchur>(
      new SetUpSchurImpl(m_field, sub_dm, is_sub, ao_up));
}