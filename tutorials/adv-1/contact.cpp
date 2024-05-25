/**
 * \file contact.cpp
 * \CONTACT contact.cpp
 *
 * CONTACT of contact problem
 *
 * @copyright Copyright (c) 2023
 */

/* The above code is a preprocessor directive in C++ that checks if the macro
"EXECUTABLE_DIMENSION" has been defined. If it has not been defined, it is set
to 3" */
#ifndef EXECUTABLE_DIMENSION
#define EXECUTABLE_DIMENSION 3
#endif

#ifndef SCHUR_ASSEMBLE
#define SCHUR_ASSEMBLE 0
#endif

#include <MoFEM.hpp>
#include <MatrixFunction.hpp>

#ifdef PYTHON_SDF
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/numpy.hpp>
namespace bp = boost::python;
namespace np = boost::python::numpy;
#endif

using namespace MoFEM;

constexpr AssemblyType AT =
    (SCHUR_ASSEMBLE) ? AssemblyType::BLOCK_SCHUR
                     : AssemblyType::PETSC; //< selected assembly type
constexpr IntegrationType IT =
    IntegrationType::GAUSS; //< selected integration type

template <int DIM> struct ElementsAndOps;

template <> struct ElementsAndOps<2> : PipelineManager::ElementsAndOpsByDim<2> {
  static constexpr FieldSpace CONTACT_SPACE = HCURL;
};

template <> struct ElementsAndOps<3> : PipelineManager::ElementsAndOpsByDim<3> {
  static constexpr FieldSpace CONTACT_SPACE = HDIV;
};

constexpr FieldSpace ElementsAndOps<2>::CONTACT_SPACE;
constexpr FieldSpace ElementsAndOps<3>::CONTACT_SPACE;

constexpr int SPACE_DIM =
    EXECUTABLE_DIMENSION; //< Space dimension of problem, mesh

/* The above code is defining an alias `EntData` for the type
`EntitiesFieldData::EntData`. This is a C++ syntax for creating a new name for
an existing type. */
using EntData = EntitiesFieldData::EntData;
using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainEleOp = DomainEle::UserDataOperator;
using BoundaryEle = ElementsAndOps<SPACE_DIM>::BoundaryEle;
using BoundaryEleOp = BoundaryEle::UserDataOperator;

//! [Specialisation for assembly]

constexpr FieldSpace CONTACT_SPACE = ElementsAndOps<SPACE_DIM>::CONTACT_SPACE;

//! [Operators used for contact]
using OpSpringLhs = FormsIntegrators<BoundaryEleOp>::Assembly<AT>::BiLinearForm<
    IT>::OpMass<1, SPACE_DIM>;
using OpSpringRhs = FormsIntegrators<BoundaryEleOp>::Assembly<AT>::LinearForm<
    IT>::OpBaseTimesVector<1, SPACE_DIM, 1>;
//! [Operators used for contact]

PetscBool is_quasi_static = PETSC_TRUE;

int order = 2;
int contact_order = 2;
int geom_order = 1;
double young_modulus = 100;
double poisson_ratio = 0.25;
double rho = 0.0;
double spring_stiffness = 0.0;
double vis_spring_stiffness = 0.0;
double alpha_damping = 0;

double scale = 1.;

PetscBool is_axisymmetric = PETSC_FALSE; //< Axisymmetric model

// ##define HENCKY_SMALL_STRAIN

int atom_test = 0;

namespace ContactOps {
double cn_contact = 0.1;
}; // namespace ContactOps

#include <HenckyOps.hpp>
using namespace HenckyOps;
#include <ContactOps.hpp>
#include <PostProcContact.hpp>
#include <ContactNaturalBC.hpp>

#ifdef WITH_MODULE_MFRONT_INTERFACE
#include <MFrontMoFEMInterface.hpp>
#endif

using DomainRhsBCs = NaturalBC<DomainEleOp>::Assembly<AT>::LinearForm<IT>;
using OpDomainRhsBCs =
    DomainRhsBCs::OpFlux<ContactOps::DomainBCs, 1, SPACE_DIM>;
using BoundaryRhsBCs = NaturalBC<BoundaryEleOp>::Assembly<AT>::LinearForm<IT>;
using OpBoundaryRhsBCs =
    BoundaryRhsBCs::OpFlux<ContactOps::BoundaryBCs, 1, SPACE_DIM>;
using BoundaryLhsBCs = NaturalBC<BoundaryEleOp>::Assembly<AT>::BiLinearForm<IT>;
using OpBoundaryLhsBCs =
    BoundaryLhsBCs::OpFlux<ContactOps::BoundaryBCs, 1, SPACE_DIM>;

using namespace ContactOps;

struct Contact {

  Contact(MoFEM::Interface &m_field) : mField(m_field) {}

  MoFEMErrorCode runProblem();

private:
  MoFEM::Interface &mField;

  MoFEMErrorCode setupProblem();
  MoFEMErrorCode createCommonData();
  MoFEMErrorCode bC();
  MoFEMErrorCode OPs();
  MoFEMErrorCode tsSolve();
  MoFEMErrorCode checkResults();

  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uXScatter;
  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uYScatter;
  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uZScatter;

  boost::shared_ptr<GenericElementInterface> mfrontInterface;
  boost::shared_ptr<Monitor> monitorPtr;

#ifdef PYTHON_SDF
  boost::shared_ptr<SDFPython> sdfPythonPtr;
#endif

  struct ScaledTimeScale : public MoFEM::TimeScale {
    using MoFEM::TimeScale::TimeScale;
    double getScale(const double time) {
      return scale * MoFEM::TimeScale::getScale(time);
    };
  };
};

//! [Run problem]
MoFEMErrorCode Contact::runProblem() {
  MoFEMFunctionBegin;
  CHKERR setupProblem();
  CHKERR createCommonData();
  CHKERR bC();
  CHKERR OPs();
  CHKERR tsSolve();
  CHKERR checkResults();
  MoFEMFunctionReturn(0);
}
//! [Run problem]

//! [Set up problem]
MoFEMErrorCode Contact::setupProblem() {
  MoFEMFunctionBegin;
  Simple *simple = mField.getInterface<Simple>();

  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-contact_order", &contact_order,
                            PETSC_NULL);
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-geom_order", &geom_order,
                            PETSC_NULL);

  MOFEM_LOG("CONTACT", Sev::inform) << "Order " << order;
  if (contact_order != order)
    MOFEM_LOG("CONTACT", Sev::inform) << "Contact order " << contact_order;
  MOFEM_LOG("CONTACT", Sev::inform) << "Geom order " << geom_order;

  // Select base
  enum bases { AINSWORTH, DEMKOWICZ, LASBASETOPT };
  const char *list_bases[LASBASETOPT] = {"ainsworth", "demkowicz"};
  PetscInt choice_base_value = AINSWORTH;
  CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                              LASBASETOPT, &choice_base_value, PETSC_NULL);

  FieldApproximationBase base;
  switch (choice_base_value) {
  case AINSWORTH:
    base = AINSWORTH_LEGENDRE_BASE;
    MOFEM_LOG("CONTACT", Sev::inform)
        << "Set AINSWORTH_LEGENDRE_BASE for displacements";
    break;
  case DEMKOWICZ:
    base = DEMKOWICZ_JACOBI_BASE;
    MOFEM_LOG("CONTACT", Sev::inform)
        << "Set DEMKOWICZ_JACOBI_BASE for displacements";
    break;
  default:
    base = LASTBASE;
    break;
  }

  // Note: For tets we have only H1 Ainsworth base, for Hex we have only H1
  // Demkowicz base. We need to implement Demkowicz H1 base on tet.
  CHKERR simple->addDomainField("U", H1, base, SPACE_DIM);
  CHKERR simple->addBoundaryField("U", H1, base, SPACE_DIM);
  CHKERR simple->addDomainField("SIGMA", CONTACT_SPACE, DEMKOWICZ_JACOBI_BASE,
                                SPACE_DIM);
  CHKERR simple->addBoundaryField("SIGMA", CONTACT_SPACE, DEMKOWICZ_JACOBI_BASE,
                                  SPACE_DIM);
  CHKERR simple->addDataField("GEOMETRY", H1, base, SPACE_DIM);

  CHKERR simple->setFieldOrder("U", order);
  CHKERR simple->setFieldOrder("GEOMETRY", geom_order);

  CHKERR simple->addDataField("SP", H1, base, SPACE_DIM);
  CHKERR simple->setFieldOrder("SP", order);

  auto get_skin = [&]() {
    Range body_ents;
    CHKERR mField.get_moab().get_entities_by_dimension(0, SPACE_DIM, body_ents);
    Skinner skin(&mField.get_moab());
    Range skin_ents;
    CHKERR skin.find_skin(0, body_ents, false, skin_ents);
    return skin_ents;
  };

  auto filter_blocks = [&](auto skin) {
    bool is_contact_block = false;
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
  int sigma_order = std::max(order, contact_order) - 1;
  CHKERR simple->setFieldOrder("SIGMA", sigma_order, &boundary_ents);

  if (contact_order > order) {
    Range ho_ents;
    if constexpr (SPACE_DIM == 3) {
      CHKERR mField.get_moab().get_adjacencies(boundary_ents, 1, false, ho_ents,
                                               moab::Interface::UNION);
    } else {
      ho_ents = boundary_ents;
    }
    CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(ho_ents);
    CHKERR simple->setFieldOrder("U", contact_order, &ho_ents);
    CHKERR mField.getInterface<CommInterface>()->synchroniseFieldEntities("U");
  }

  CHKERR simple->setUp();

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
} //! [Set up problem]

//! [Create common data]
MoFEMErrorCode Contact::createCommonData() {
  MoFEMFunctionBegin;

  PetscBool use_mfront = PETSC_FALSE;
  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-use_mfront", &use_mfront,
                             PETSC_NULL);
  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-is_axisymmetric",
                             &is_axisymmetric, PETSC_NULL);
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-atom_test", &atom_test,
                            PETSC_NULL);

  CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-scale", &scale, PETSC_NULL);
  CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-young_modulus", &young_modulus,
                               PETSC_NULL);
  CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-poisson_ratio", &poisson_ratio,
                               PETSC_NULL);
  CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-rho", &rho, PETSC_NULL);
  CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-cn", &cn_contact, PETSC_NULL);
  CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-spring_stiffness",
                               &spring_stiffness, PETSC_NULL);
  CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-vis_spring_stiffness",
                               &vis_spring_stiffness, PETSC_NULL);
  CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-alpha_damping", &alpha_damping,
                               PETSC_NULL);

  if (!mfrontInterface) {
    MOFEM_LOG("CONTACT", Sev::inform) << "Young modulus " << young_modulus;
    MOFEM_LOG("CONTACT", Sev::inform) << "Poisson_ratio " << poisson_ratio;
  } else {
    MOFEM_LOG("CONTACT", Sev::inform) << "Using MFront for material model";
  }

  MOFEM_LOG("CONTACT", Sev::inform) << "Density " << rho;
  MOFEM_LOG("CONTACT", Sev::inform) << "cn_contact " << cn_contact;
  MOFEM_LOG("CONTACT", Sev::inform) << "Spring stiffness " << spring_stiffness;
  MOFEM_LOG("CONTACT", Sev::inform)
      << "Vis spring_stiffness " << vis_spring_stiffness;

  MOFEM_LOG("CONTACT", Sev::inform) << "alpha_damping " << alpha_damping;

  PetscBool use_scale = PETSC_FALSE;
  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-use_scale", &use_scale,
                             PETSC_NULL);
  if (use_scale) {
    scale /= young_modulus;
  }

  MOFEM_LOG("CONTACT", Sev::inform) << "Scale " << scale;

  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-is_quasi_static",
                             &is_quasi_static, PETSC_NULL);
  MOFEM_LOG("CONTACT", Sev::inform)
      << "Is quasi-static: " << (is_quasi_static ? "true" : "false");

#ifdef PYTHON_SDF
  char sdf_file_name[255];
  CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-sdf_file",
                               sdf_file_name, 255, PETSC_NULL);

  sdfPythonPtr = boost::make_shared<SDFPython>();
  CHKERR sdfPythonPtr->sdfInit(sdf_file_name);
  sdfPythonWeakPtr = sdfPythonPtr;
#endif

  if (is_axisymmetric) {
    if (SPACE_DIM == 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Use executable contact_2d with axisymmetric model");
    } else {
      if (!use_mfront) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                "Axisymmetric model is only available with MFront (set "
                "use_mfront to 1)");
      } else {
        MOFEM_LOG("CONTACT", Sev::inform) << "Using axisymmetric model";
      }
    }
  } else {
    if (SPACE_DIM == 2) {
      MOFEM_LOG("CONTACT", Sev::inform) << "Using plane strain model";
    }
  }

  if (use_mfront) {
#ifndef WITH_MODULE_MFRONT_INTERFACE
    SETERRQ(
        PETSC_COMM_SELF, MOFEM_NOT_FOUND,
        "MFrontInterface module was not found while use_mfront was set to 1");
#else
    if (SPACE_DIM == 3) {
      mfrontInterface =
          boost::make_shared<MFrontMoFEMInterface<TRIDIMENSIONAL>>(
              mField, "U", "GEOMETRY", true, is_quasi_static);
    } else if (SPACE_DIM == 2) {
      if (is_axisymmetric) {
        mfrontInterface =
            boost::make_shared<MFrontMoFEMInterface<AXISYMMETRICAL>>(
                mField, "U", "GEOMETRY", true, is_quasi_static);
      } else {
        mfrontInterface = boost::make_shared<MFrontMoFEMInterface<PLANESTRAIN>>(
            mField, "U", "GEOMETRY", true, is_quasi_static);
      }
    }
#endif
    CHKERR mfrontInterface->getCommandLineParameters();
  }

  Simple *simple = mField.getInterface<Simple>();
  auto dm = simple->getDM();
  monitorPtr =
      boost::make_shared<Monitor>(dm, scale, mfrontInterface, is_axisymmetric);

  if (use_mfront) {
    mfrontInterface->setMonitorPtr(monitorPtr);
  }

  MoFEMFunctionReturn(0);
}
//! [Create common data]

//! [Boundary condition]
MoFEMErrorCode Contact::bC() {
  MoFEMFunctionBegin;
  auto bc_mng = mField.getInterface<BcManager>();
  auto simple = mField.getInterface<Simple>();

  for (auto f : {"U", "SIGMA"}) {
    CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(),
                                             "REMOVE_X", f, 0, 0);
    CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(),
                                             "REMOVE_Y", f, 1, 1);
    CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(),
                                             "REMOVE_Z", f, 2, 2);
    CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(),
                                             "REMOVE_ALL", f, 0, 3);
  }

  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "FIX_X",
                                           "SIGMA", 0, 0, false, true);
  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "FIX_Y",
                                           "SIGMA", 1, 1, false, true);
  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "FIX_Z",
                                           "SIGMA", 2, 2, false, true);
  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "FIX_ALL",
                                           "SIGMA", 0, 3, false, true);
  CHKERR bc_mng->removeBlockDOFsOnEntities(
      simple->getProblemName(), "NO_CONTACT", "SIGMA", 0, 3, false, true);

  // Note remove has to be always before push. Then node marking will be
  // corrupted.
  CHKERR bc_mng->pushMarkDOFsOnEntities<DisplacementCubitBcData>(
      simple->getProblemName(), "U");

  MoFEMFunctionReturn(0);
}
//! [Boundary condition]

//! [Push operators to pip]
MoFEMErrorCode Contact::OPs() {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  auto *pip_mng = mField.getInterface<PipelineManager>();
  auto bc_mng = mField.getInterface<BcManager>();
  auto time_scale = boost::make_shared<ScaledTimeScale>();
  auto body_force_time_scale =
      boost::make_shared<ScaledTimeScale>("body_force_hist.txt");

  auto integration_rule_vol = [](int, int, int approx_order) {
    return 2 * approx_order + geom_order - 1;
  };
  auto integration_rule_boundary = [](int, int, int approx_order) {
    return 2 * approx_order + geom_order - 1;
  };

  auto add_domain_base_ops = [&](auto &pip) {
    MoFEMFunctionBegin;
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1, HDIV},
                                                          "GEOMETRY");
    MoFEMFunctionReturn(0);
  };

  auto henky_common_data_ptr = boost::make_shared<HenckyOps::CommonData>();
  henky_common_data_ptr->matDPtr = boost::make_shared<MatrixDouble>();
  henky_common_data_ptr->matGradPtr = boost::make_shared<MatrixDouble>();

  auto add_domain_ops_lhs = [&](auto &pip) {
    MoFEMFunctionBegin;

    //! [Only used for dynamics]
    using OpMass = FormsIntegrators<DomainEleOp>::Assembly<AT>::BiLinearForm<
        GAUSS>::OpMass<1, SPACE_DIM>;
    //! [Only used for dynamics]
    if (is_quasi_static == PETSC_FALSE) {

      auto *pip_mng = mField.getInterface<PipelineManager>();
      auto fe_domain_lhs = pip_mng->getDomainLhsFE();

      auto get_inertia_and_mass_damping =
          [this, fe_domain_lhs](const double, const double, const double) {
            return (rho * scale) * fe_domain_lhs->ts_aa +
                   (alpha_damping * scale) * fe_domain_lhs->ts_a;
          };
      pip.push_back(new OpMass("U", "U", get_inertia_and_mass_damping));
    } else {

      auto *pip_mng = mField.getInterface<PipelineManager>();
      auto fe_domain_lhs = pip_mng->getDomainLhsFE();

      auto get_inertia_and_mass_damping =
          [this, fe_domain_lhs](const double, const double, const double) {
            return (alpha_damping * scale) * fe_domain_lhs->ts_a;
          };
      pip.push_back(new OpMass("U", "U", get_inertia_and_mass_damping));
    }

    if (!mfrontInterface) {
      CHKERR HenckyOps::opFactoryDomainLhs<SPACE_DIM, AT, IT, DomainEleOp>(
          mField, pip, "U", "MAT_ELASTIC", Sev::verbose, scale);
    } else {
      CHKERR mfrontInterface->opFactoryDomainLhs(pip);
    }

    MoFEMFunctionReturn(0);
  };

  auto add_domain_ops_rhs = [&](auto &pip) {
    MoFEMFunctionBegin;

    CHKERR DomainRhsBCs::AddFluxToPipeline<OpDomainRhsBCs>::add(
        pip, mField, "U", {body_force_time_scale}, Sev::inform);

    //! [Only used for dynamics]
    using OpInertiaForce = FormsIntegrators<DomainEleOp>::Assembly<
        AT>::LinearForm<IT>::OpBaseTimesVector<1, SPACE_DIM, 1>;
    //! [Only used for dynamics]

    // only in case of dynamics
    if (is_quasi_static == PETSC_FALSE) {
      auto mat_acceleration = boost::make_shared<MatrixDouble>();
      pip.push_back(new OpCalculateVectorFieldValuesDotDot<SPACE_DIM>(
          "U", mat_acceleration));
      pip.push_back(
          new OpInertiaForce("U", mat_acceleration, [](double, double, double) {
            return rho * scale;
          }));
    }

    // only in case of viscosity
    if (alpha_damping > 0) {
      auto mat_velocity = boost::make_shared<MatrixDouble>();
      pip.push_back(
          new OpCalculateVectorFieldValuesDot<SPACE_DIM>("U", mat_velocity));
      pip.push_back(
          new OpInertiaForce("U", mat_velocity, [](double, double, double) {
            return alpha_damping * scale;
          }));
    }

    if (!mfrontInterface) {
      CHKERR HenckyOps::opFactoryDomainRhs<SPACE_DIM, AT, IT, DomainEleOp>(
          mField, pip, "U", "MAT_ELASTIC", Sev::inform, scale);
    } else {
      CHKERR mfrontInterface->opFactoryDomainRhs(pip);
    }

    CHKERR ContactOps::opFactoryDomainRhs<SPACE_DIM, AT, IT, DomainEleOp>(
        pip, "SIGMA", "U", is_axisymmetric);

    MoFEMFunctionReturn(0);
  };

  auto add_boundary_base_ops = [&](auto &pip) {
    MoFEMFunctionBegin;
    CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(pip, {HDIV},
                                                              "GEOMETRY");
    // We have to integrate on curved face geometry, thus integration weight
    // have to adjusted.
    pip.push_back(new OpSetHOWeightsOnSubDim<SPACE_DIM>());
    MoFEMFunctionReturn(0);
  };

  auto add_boundary_ops_lhs = [&](auto &pip) {
    MoFEMFunctionBegin;

    //! [Operators used for contact]
    using OpSpringLhs = FormsIntegrators<BoundaryEleOp>::Assembly<
        AT>::BiLinearForm<IT>::OpMass<1, SPACE_DIM>;
    //! [Operators used for contact]

    // Add Natural BCs to LHS
    CHKERR BoundaryLhsBCs::AddFluxToPipeline<OpBoundaryLhsBCs>::add(
        pip, mField, "U", Sev::inform);

    if (spring_stiffness > 0 || vis_spring_stiffness > 0) {

      auto *pip_mng = mField.getInterface<PipelineManager>();
      auto fe_boundary_lhs = pip_mng->getBoundaryLhsFE();

      pip.push_back(new OpSpringLhs(
          "U", "U",

          [this, fe_boundary_lhs](double, double, double) {
            return spring_stiffness * scale +
                   (vis_spring_stiffness * scale) * fe_boundary_lhs->ts_a;
          }

          ));
    }

    CHKERR
    ContactOps::opFactoryBoundaryLhs<SPACE_DIM, AT, GAUSS, BoundaryEleOp>(
        pip, "SIGMA", "U", is_axisymmetric);
    CHKERR ContactOps::opFactoryBoundaryToDomainLhs<SPACE_DIM, AT, GAUSS,
                                                    DomainEle>(
        mField, pip, simple->getDomainFEName(), "SIGMA", "U", "GEOMETRY",
        integration_rule_vol, is_axisymmetric);

    MoFEMFunctionReturn(0);
  };

  auto add_boundary_ops_rhs = [&](auto &pip) {
    MoFEMFunctionBegin;

    //! [Operators used for contact]
    using OpSpringRhs = FormsIntegrators<BoundaryEleOp>::Assembly<
        AT>::LinearForm<IT>::OpBaseTimesVector<1, SPACE_DIM, 1>;
    //! [Operators used for contact]

    // Add Natural BCs to RHS
    CHKERR BoundaryRhsBCs::AddFluxToPipeline<OpBoundaryRhsBCs>::add(
        pip, mField, "U", {time_scale}, Sev::inform);

    if (spring_stiffness > 0 || vis_spring_stiffness > 0) {
      auto u_disp = boost::make_shared<MatrixDouble>();
      auto dot_u_disp = boost::make_shared<MatrixDouble>();
      pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_disp));
      pip.push_back(
          new OpCalculateVectorFieldValuesDot<SPACE_DIM>("U", dot_u_disp));
      pip.push_back(
          new OpSpringRhs("U", u_disp, [this](double, double, double) {
            return spring_stiffness * scale;
          }));
      pip.push_back(
          new OpSpringRhs("U", dot_u_disp, [this](double, double, double) {
            return vis_spring_stiffness * scale;
          }));
    }

    CHKERR
    ContactOps::opFactoryBoundaryRhs<SPACE_DIM, AT, GAUSS, BoundaryEleOp>(
        pip, "SIGMA", "U", is_axisymmetric);

    MoFEMFunctionReturn(0);
  };

  CHKERR add_domain_base_ops(pip_mng->getOpDomainLhsPipeline());
  CHKERR add_domain_base_ops(pip_mng->getOpDomainRhsPipeline());
  CHKERR add_domain_ops_lhs(pip_mng->getOpDomainLhsPipeline());
  CHKERR add_domain_ops_rhs(pip_mng->getOpDomainRhsPipeline());

  CHKERR add_boundary_base_ops(pip_mng->getOpBoundaryLhsPipeline());
  CHKERR add_boundary_base_ops(pip_mng->getOpBoundaryRhsPipeline());
  CHKERR add_boundary_ops_lhs(pip_mng->getOpBoundaryLhsPipeline());
  CHKERR add_boundary_ops_rhs(pip_mng->getOpBoundaryRhsPipeline());

  if (mfrontInterface) {
    CHKERR mfrontInterface->setUpdateElementVariablesOperators();
  }

  CHKERR pip_mng->setDomainRhsIntegrationRule(integration_rule_vol);
  CHKERR pip_mng->setDomainLhsIntegrationRule(integration_rule_vol);
  CHKERR pip_mng->setBoundaryRhsIntegrationRule(integration_rule_boundary);
  CHKERR pip_mng->setBoundaryLhsIntegrationRule(integration_rule_boundary);

  MoFEMFunctionReturn(0);
}
//! [Push operators to pip]

//! [Solve]
struct SetUpSchur {
  static boost::shared_ptr<SetUpSchur>
  createSetUpSchur(MoFEM::Interface &m_field);
  virtual MoFEMErrorCode setUp(SmartPetscObj<TS> solver) = 0;

protected:
  SetUpSchur() = default;
};

MoFEMErrorCode Contact::tsSolve() {
  MoFEMFunctionBegin;

  Simple *simple = mField.getInterface<Simple>();
  PipelineManager *pip_mng = mField.getInterface<PipelineManager>();
  ISManager *is_manager = mField.getInterface<ISManager>();

  auto set_section_monitor = [&](auto solver) {
    MoFEMFunctionBegin;
    SNES snes;
    CHKERR TSGetSNES(solver, &snes);
    PetscViewerAndFormat *vf;
    CHKERR PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,
                                      PETSC_VIEWER_DEFAULT, &vf);
    CHKERR SNESMonitorSet(
        snes,
        (MoFEMErrorCode(*)(SNES, PetscInt, PetscReal, void *))SNESMonitorFields,
        vf, (MoFEMErrorCode(*)(void **))PetscViewerAndFormatDestroy);
    MoFEMFunctionReturn(0);
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
    monitorPtr->setScatterVectors(uXScatter, uYScatter, uZScatter);
    boost::shared_ptr<ForcesAndSourcesCore> null;
    CHKERR DMMoFEMTSSetMonitor(dm, solver, simple->getDomainFEName(),
                               monitorPtr, null, null);
    MoFEMFunctionReturn(0);
  };

  auto set_essential_bc = [&]() {
    MoFEMFunctionBegin;
    // This is low level pushing finite elements (pipelines) to solver
    auto ts_ctx_ptr = getDMTsCtx(simple->getDM());
    auto pre_proc_ptr = boost::make_shared<FEMethod>();
    auto post_proc_rhs_ptr = boost::make_shared<FEMethod>();
    auto post_proc_lhs_ptr = boost::make_shared<FEMethod>();

    // Add boundary condition scaling
    auto time_scale = boost::make_shared<TimeScale>();

    auto get_bc_hook_rhs = [&]() {
      EssentialPreProc<DisplacementCubitBcData> hook(mField, pre_proc_ptr,
                                                     {time_scale}, false);
      return hook;
    };
    pre_proc_ptr->preProcessHook = get_bc_hook_rhs();

    auto get_post_proc_hook_rhs = [&]() {
      return EssentialPostProcRhs<DisplacementCubitBcData>(
          mField, post_proc_rhs_ptr, 1.);
    };
    auto get_post_proc_hook_lhs = [&]() {
      return EssentialPostProcLhs<DisplacementCubitBcData>(
          mField, post_proc_lhs_ptr, 1.);
    };
    post_proc_rhs_ptr->postProcessHook = get_post_proc_hook_rhs();

    ts_ctx_ptr->getPreProcessIFunction().push_front(pre_proc_ptr);
    ts_ctx_ptr->getPreProcessIJacobian().push_front(pre_proc_ptr);
    ts_ctx_ptr->getPostProcessIFunction().push_back(post_proc_rhs_ptr);
    post_proc_lhs_ptr->postProcessHook = get_post_proc_hook_lhs();
    ts_ctx_ptr->getPostProcessIJacobian().push_back(post_proc_lhs_ptr);
    MoFEMFunctionReturn(0);
  };

  // Set up Schur preconditioner
  auto set_schur_pc = [&](auto solver) {
    boost::shared_ptr<SetUpSchur> schur_ptr;
    if (AT == AssemblyType::BLOCK_SCHUR) {
      // Set up Schur preconditioner
      schur_ptr = SetUpSchur::createSetUpSchur(mField);
      CHK_MOAB_THROW(schur_ptr->setUp(solver), "SetUpSchur::setUp");
    }
    return schur_ptr;
  };

  auto dm = simple->getDM();
  auto D = createDMVector(dm);

  ContactOps::CommonData::createTotalTraction(mField);

  uXScatter = scatter_create(D, 0);
  uYScatter = scatter_create(D, 1);
  if (SPACE_DIM == 3)
    uZScatter = scatter_create(D, 2);

  // Add extra finite elements to SNES solver pipelines to resolve essential
  // boundary conditions
  CHKERR set_essential_bc();

  if (is_quasi_static == PETSC_TRUE) {
    auto solver = pip_mng->createTSIM();
    CHKERR TSSetFromOptions(solver);
    auto schur_pc_ptr = set_schur_pc(solver);

    auto D = createDMVector(dm);
    CHKERR set_section_monitor(solver);
    CHKERR set_time_monitor(dm, solver);
    CHKERR TSSetSolution(solver, D);
    CHKERR TSSetUp(solver);
    CHKERR TSSolve(solver, NULL);
  } else {
    auto solver = pip_mng->createTSIM2();
    CHKERR TSSetFromOptions(solver);
    auto schur_pc_ptr = set_schur_pc(solver);

    auto dm = simple->getDM();
    auto D = createDMVector(dm);
    auto DD = vectorDuplicate(D);

    CHKERR set_section_monitor(solver);
    CHKERR set_time_monitor(dm, solver);
    CHKERR TS2SetSolution(solver, D, DD);
    CHKERR TSSetUp(solver);
    CHKERR TSSolve(solver, NULL);
  }

  MoFEMFunctionReturn(0);
}
//! [Solve]

//! [Check]
MoFEMErrorCode Contact::checkResults() {
  MoFEMFunctionBegin;
  if (atom_test && !mField.get_comm_rank()) {
    const double *t_ptr;
    CHKERR VecGetArrayRead(ContactOps::CommonData::totalTraction, &t_ptr);
    double hertz_force;
    double fem_force;
    double tol = 1e-3;
    switch (atom_test) {
    case 1: // plane stress
      hertz_force = 3.927;
      fem_force = t_ptr[1];
      break;
    case 2: // plane strain
      hertz_force = 4.675;
      fem_force = t_ptr[1];
      break;
    case 3: // 3D
      hertz_force = 3.968;
      fem_force = t_ptr[2];
    case 4: // axisymmetric
      tol = 5e3;
    case 5: // axisymmetric
      hertz_force = 15.873;
      fem_force = t_ptr[1];
      break;
    case 6: // wavy 2d
      hertz_force = 0.374;
      fem_force = t_ptr[1];
      break;
    case 7: // wavy 3d
      hertz_force = 0.5289;
      fem_force = t_ptr[2];
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
               "atom test %d does not exist", atom_test);
    }
    if (fabs(fem_force - hertz_force) / hertz_force > tol) {
      SETERRQ3(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
               "atom test %d diverged! %3.4e != %3.4e", atom_test, fem_force,
               hertz_force);
    }
    CHKERR VecRestoreArrayRead(ContactOps::CommonData::totalTraction, &t_ptr);
  }

  ContactOps::CommonData::totalTraction.reset();

  MoFEMFunctionReturn(0);
}
//! [Check]

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

#ifdef PYTHON_SDF
  Py_Initialize();
  np::initialize();
#endif

  // Initialisation of MoFEM/PETSc and MOAB data structures
  const char param_file[] = "param_file.petsc";
  MoFEM::Core::Initialize(&argc, &argv, param_file, help);

  // Add logging channel for CONTACT
  auto core_log = logging::core::get();
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmWorld(), "CONTACT"));
  LogManager::setLog("CONTACT");
  MOFEM_LOG_TAG("CONTACT", "Indent");

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
    CHKERR simple->loadFile("");
    //! [Load mesh]

    //! [CONTACT]
    Contact ex(m_field);
    CHKERR ex.runProblem();
    //! [CONTACT]
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

#ifdef PYTHON_SDF
  if (Py_FinalizeEx() < 0) {
    exit(120);
  }
#endif

  return 0;
}

struct SetUpSchurImpl : public SetUpSchur {

  SetUpSchurImpl(MoFEM::Interface &m_field) : SetUpSchur(), mField(m_field) {}

  virtual ~SetUpSchurImpl() {}

  MoFEMErrorCode setUp(SmartPetscObj<TS> solver);

private:
  MoFEMErrorCode createSubDM();
  MoFEMErrorCode setOperator();
  MoFEMErrorCode setPC(PC pc);
  MoFEMErrorCode setDiagonalPC(PC pc);

  MoFEM::Interface &mField;

  SmartPetscObj<Mat> S;
  SmartPetscObj<DM> schurDM;
  SmartPetscObj<DM> blockDM;
};

MoFEMErrorCode SetUpSchurImpl::setUp(SmartPetscObj<TS> solver) {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  auto pip = mField.getInterface<PipelineManager>();

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

    MOFEM_LOG("CONTACT", Sev::inform) << "Setup Schur pc";

    if (S) {
      CHK_THROW_MESSAGE(
          MOFEM_DATA_INCONSISTENCY,
          "It is expected that Schur matrix is not allocated. This is "
          "possible only if PC is set up twice");
    }

    CHKERR createSubDM();

    // Add data to DM storage
    S = createDMMatrix(schurDM);
    CHKERR MatSetBlockSize(S, SPACE_DIM);
    // CHKERR MatSetOption(S, MAT_SYMMETRIC, PETSC_TRUE);

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

    CHKERR setOperator();
    CHKERR setPC(pc);
    CHKERR TSSetUp(solver);
    CHKERR KSPSetUp(ksp);
    CHKERR setDiagonalPC(pc);

  } else {
    MOFEM_LOG("CONTACT", Sev::inform) << "No Schur PC";
    pip->getOpBoundaryLhsPipeline().push_front(createOpSchurAssembleBegin());
    pip->getOpBoundaryLhsPipeline().push_back(
        createOpSchurAssembleEnd({}, {}, {}, {}, {}, false));
    pip->getOpDomainLhsPipeline().push_front(createOpSchurAssembleBegin());
    pip->getOpDomainLhsPipeline().push_back(
        createOpSchurAssembleEnd({}, {}, {}, {}, {}, false));
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SetUpSchurImpl::createSubDM() {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();

  auto create_dm = [&](const char *name, const char *field_name) {
    auto dm = createDM(mField.get_comm(), "DMMOFEM");
    auto create_dm_imp = [&]() {
      MoFEMFunctionBegin;
      CHKERR DMMoFEMCreateSubDM(dm, simple->getDM(), name);
      CHKERR DMMoFEMSetSquareProblem(dm, PETSC_TRUE);
      CHKERR DMMoFEMAddElement(dm, simple->getDomainFEName());
      CHKERR DMMoFEMAddSubFieldRow(dm, field_name);
      CHKERR DMMoFEMAddSubFieldCol(dm, field_name);
      CHKERR DMSetUp(dm);
      MoFEMFunctionReturn(0);
    };
    CHK_THROW_MESSAGE(
        create_dm_imp(),
        "Error in creating schurDM. It is possible that schurDM is "
        "already created");
    return dm;
  };

  // Note: here we can make block with bubbles of "U" and "SIGMA" fields. See
  // vec-0 where bubbles are added.

  schurDM = create_dm("SCHUR", "U");
  blockDM = create_dm("BLOCK", "SIGMA");

  if constexpr (AT == AssemblyType::BLOCK_SCHUR) {

    auto get_nested_mat_data = [&](auto schur_dm, auto block_dm) {
      auto block_mat_data = createBlockMatStructure(
          simple->getDM(),

          {{

              simple->getDomainFEName(),

              {

                  {"U", "U"}, {"SIGMA", "U"}, {"U", "SIGMA"}, {"SIGMA", "SIGMA"}

              }}}

      );

      return getNestSchurData(

          {schurDM, blockDM}, block_mat_data,

          {"SIGMA"}, {nullptr}, true

      );
    };

    auto nested_mat_data = get_nested_mat_data(schurDM, blockDM);
    CHKERR DMMoFEMSetNestSchurData(simple->getDM(), nested_mat_data);

  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Only BLOCK_SCHUR is implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SetUpSchurImpl::setOperator() {
  MoFEMFunctionBegin;

  double eps_stab = 1e-4;
  CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-eps_stab", &eps_stab,
                               PETSC_NULL);

  using B = FormsIntegrators<BoundaryEleOp>::Assembly<
      BLOCK_PRECONDITIONER_SCHUR>::BiLinearForm<IT>;
  using OpMassStab = B::OpMass<3, SPACE_DIM * SPACE_DIM>;

  auto simple = mField.getInterface<Simple>();
  auto pip = mField.getInterface<PipelineManager>();

  // block data structure
  boost::shared_ptr<BlockStructure> block_data;
  CHKERR DMMoFEMGetBlocMatData(simple->getDM(), block_data);

  // Boundary
  auto dm_is = getDMSubData(schurDM)->getSmartRowIs();
  auto ao_up = createAOMappingIS(dm_is, PETSC_NULL);

  pip->getOpBoundaryLhsPipeline().push_front(createOpSchurAssembleBegin());
  pip->getOpBoundaryLhsPipeline().push_back(
      new OpMassStab("SIGMA", "SIGMA",
                     [eps_stab](double, double, double) { return eps_stab; }));
  pip->getOpBoundaryLhsPipeline().push_back(

      createOpSchurAssembleEnd({"SIGMA"}, {nullptr}, {ao_up}, {S}, {false},
                               false, block_data)

  );

  // Domain
  pip->getOpDomainLhsPipeline().push_front(createOpSchurAssembleBegin());
  pip->getOpDomainLhsPipeline().push_back(

      createOpSchurAssembleEnd({"SIGMA"}, {nullptr}, {ao_up}, {S}, {false},
                               false, block_data)

  );

  auto pre_proc_schur_lhs_ptr = boost::make_shared<FEMethod>();
  auto post_proc_schur_lhs_ptr = boost::make_shared<FEMethod>();

  pre_proc_schur_lhs_ptr->preProcessHook = [this]() {
    MoFEMFunctionBegin;
    CHKERR MatZeroEntries(S);
    MOFEM_LOG("CONTACT", Sev::verbose) << "Lhs Assemble Begin";
    MoFEMFunctionReturn(0);
  };

  post_proc_schur_lhs_ptr->postProcessHook = [this, ao_up,
                                              post_proc_schur_lhs_ptr]() {
    MoFEMFunctionBegin;
    MOFEM_LOG("CONTACT", Sev::verbose) << "Lhs Assemble End";
    auto print_mat_norm = [this](auto a, std::string prefix) {
      MoFEMFunctionBegin;
      double nrm;
      CHKERR MatNorm(a, NORM_FROBENIUS, &nrm);
      MOFEM_LOG("CONTACT", Sev::noisy) << prefix << " norm = " << nrm;
      MoFEMFunctionReturn(0);
    };
    CHKERR MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
    CHKERR MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
    CHKERR EssentialPostProcLhs<DisplacementCubitBcData>(
        mField, post_proc_schur_lhs_ptr, 1, S, ao_up)();
#ifndef NDEBUG
    CHKERR print_mat_norm(S, "S");
#endif // NDEBUG
    MOFEM_LOG("CONTACT", Sev::verbose) << "Lhs Assemble Finish";
    MoFEMFunctionReturn(0);
  };

  auto ts_ctx_ptr = getDMTsCtx(simple->getDM());
  ts_ctx_ptr->getPreProcessIJacobian().push_front(pre_proc_schur_lhs_ptr);
  ts_ctx_ptr->getPostProcessIJacobian().push_back(post_proc_schur_lhs_ptr);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SetUpSchurImpl::setPC(PC pc) {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  auto block_is = getDMSubData(blockDM)->getSmartRowIs();
  CHKERR PCFieldSplitSetIS(pc, NULL, block_is);
  CHKERR PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, S);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SetUpSchurImpl::setDiagonalPC(PC pc) {
  MoFEMFunctionBegin;
  KSP *subksp;
  CHKERR PCFieldSplitSchurGetSubKSP(pc, PETSC_NULL, &subksp);
  auto get_pc = [](auto ksp) {
    PC pc_raw;
    CHKERR KSPGetPC(ksp, &pc_raw);
    return SmartPetscObj<PC>(pc_raw, true); // bump reference
  };
  CHKERR setSchurMatSolvePC(get_pc(subksp[0]));
  CHKERR PetscFree(subksp);
  MoFEMFunctionReturn(0);
}

boost::shared_ptr<SetUpSchur>
SetUpSchur::createSetUpSchur(MoFEM::Interface &m_field) {
  return boost::shared_ptr<SetUpSchur>(new SetUpSchurImpl(m_field));
}
