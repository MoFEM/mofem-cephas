/** \file LoopMethods.hpp
 * \brief MoFEM interface
 *
 * Data structures for making loops over finite elements and entities in the
 * problem or MoFEM database.
 *
 */

/*
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef __LOOPMETHODS_HPP__
#define __LOOPMETHODS_HPP__

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMPetscDataMethod =
    MOFEMuuid(BitIntefaceId(PETSC_DATA_METHOD));
static const MOFEMuuid IDD_MOFEMKspMethod =
    MOFEMuuid(BitIntefaceId(KSP_METHOD));
static const MOFEMuuid IDD_MOFEMSnesMethod =
    MOFEMuuid(BitIntefaceId(SNES_METHOD));
static const MOFEMuuid IDD_MOFEMTsMethod = MOFEMuuid(BitIntefaceId(TS_METHOD));
static const MOFEMuuid IDD_MOFEMBasicMethod =
    MOFEMuuid(BitIntefaceId(BASIC_METHOD));
static const MOFEMuuid IDD_MOFEMFEMethod = MOFEMuuid(BitIntefaceId(FE_METHOD));
static const MOFEMuuid IDD_MOFEMEntityMethod =
    MOFEMuuid(BitIntefaceId(ENTITY_METHOD));
static const MOFEMuuid IDD_MOFEMDofMethod =
    MOFEMuuid(BitIntefaceId(DOF_METHOD));

struct PetscData : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  PetscData();

  virtual ~PetscData() = default;

  enum DataContext {
    CTX_SET_NONE = 0,
    CTX_SET_F = 1 << 0,
    CTX_SET_A = 1 << 1,
    CTX_SET_B = 1 << 2,
    CTX_SET_X = 1 << 3,
    CTX_SET_X_T = 1 << 4,
    CTX_SET_X_TT = 1 << 6,
    CTX_SET_TIME = 1 << 7
  };

  using Switches = std::bitset<8>;

  static constexpr Switches CtxSetNone = PetscData::Switches(CTX_SET_NONE);
  static constexpr Switches CtxSetF = PetscData::Switches(CTX_SET_F);
  static constexpr Switches CtxSetA = PetscData::Switches(CTX_SET_A);
  static constexpr Switches CtxSetB = PetscData::Switches(CTX_SET_B);
  static constexpr Switches CtxSetX = PetscData::Switches(CTX_SET_X);
  static constexpr Switches CtxSetX_T = PetscData::Switches(CTX_SET_X_T);
  static constexpr Switches CtxSetX_TT = PetscData::Switches(CTX_SET_X_TT);
  static constexpr Switches CtxSetTime = PetscData::Switches(CTX_SET_TIME);

  Switches data_ctx;

  Vec f;
  Mat A;
  Mat B;
  Vec x;
  Vec x_t;
  Vec x_tt;
};

/**
 * \brief data structure for ksp (linear solver) context
 * \ingroup mofem_loops
 *
 * Struture stores context data which are set in functions run by PETSc SNES
 * functions.
 *
 */
struct KspMethod : virtual public PetscData {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;
  /**
   * \brief pass information about context of KSP/DM for with finite element is
   * computed
   */
  enum KSPContext { CTX_SETFUNCTION, CTX_OPERATORS, CTX_KSPNONE };

  KspMethod();

  virtual ~KspMethod() = default;

  /**
   * \brief copy data form another method
   * @param  ksp ksp method
   * @return     error code
   */
  MoFEMErrorCode copyKsp(const KspMethod &ksp);

  KSPContext ksp_ctx; ///< Context
  KSP ksp;            ///< KSP solver

  Vec &ksp_f;
  Mat &ksp_A;
  Mat &ksp_B;
};

/**
 * \brief data structure for snes (nonlinear solver) context
 * \ingroup mofem_loops
 *
 * Structure stores context data which are set in functions run by PETSc SNES
 * functions.
 *
 */
struct SnesMethod : virtual protected PetscData {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  enum SNESContext { CTX_SNESSETFUNCTION, CTX_SNESSETJACOBIAN, CTX_SNESNONE };

  SnesMethod();

  virtual ~SnesMethod() = default;

  /**
   * \brief Copy snes data
   */
  MoFEMErrorCode copySnes(const SnesMethod &snes);

  SNESContext snes_ctx;

  /**
   * @deprecated Avoid using values by hand.
   */
  DEPRECATED inline MoFEMErrorCode setSnesCtx(SNESContext ctx);

  SNES snes;   ///< snes solver
  Vec &snes_x; ///< state vector
  Vec &snes_f; ///< residual
  Mat &snes_A; ///< jacobian matrix
  Mat &snes_B; ///< preconditioner of jacobian matrix
};

MoFEMErrorCode SnesMethod::setSnesCtx(SNESContext ctx) {
  snes_ctx = ctx;
  return 0;
}

/**
 * \brief data structure for TS (time stepping) context
 * \ingroup mofem_loops
 *
 * Structure stores context data which are set in functions run by PETSc Time
 * Stepping functions.
 */
struct TSMethod : virtual protected PetscData {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  enum TSContext {
    CTX_TSSETRHSFUNCTION,
    CTX_TSSETRHSJACOBIAN,
    CTX_TSSETIFUNCTION,
    CTX_TSSETIJACOBIAN,
    CTX_TSTSMONITORSET,
    CTX_TSNONE
  };

  TSMethod();

  virtual ~TSMethod() = default;

  /// \brief Copy TS solver data
  MoFEMErrorCode copyTs(const TSMethod &ts);

  TS ts; ///< time solver

  TSContext ts_ctx;

  /**
   * @deprecated Avoid using values by hand.
   */
  DEPRECATED inline MoFEMErrorCode setTsCtx(TSContext ctx);

  PetscInt ts_step; ///< time step
  PetscReal ts_a;   ///< shift for U_tt (see PETSc Time Solver)
  PetscReal ts_v;   ///< shift for U_t shift for U_t
  PetscReal ts_t;   ///< time

  Vec &ts_u;    ///< state vector
  Vec &ts_u_t;  ///< time derivative of state vector
  Vec &ts_u_tt; ///< second time derivative of state vector
  Vec &ts_F;    ///< residual vector

  Mat &ts_A; ///< Jacobian of G(U) = F(t,U,W+v*U,W'+a*U), equivalent to dF/dU +
             ///< v*dF/dU_t + a*dF/dU_tt
  Mat &ts_B; ///< Preconditioner for ts_A
};

MoFEMErrorCode TSMethod::setTsCtx(TSContext ctx) {
  ts_ctx = ctx;
  return 0;
}

/**
 * \brief Data structure to exchange data between mofem and User Loop Methods.
 * \ingroup mofem_loops
 *
 * It allows to exchange data between MoFEM and user functions. It stores
 * information about multi-indices.
 *
 */
struct BasicMethod : public KspMethod, SnesMethod, TSMethod {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const {
    if (uuid == IDD_MOFEMBasicMethod) {
      *iface = const_cast<BasicMethod *>(this);
      MoFEMFunctionReturnHot(0);
    }
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  }

  BasicMethod();
  virtual ~BasicMethod() = default;

  /**
   * @brief number currently of processed method
   */
  int nInTheLoop;

  /**
   * @brief local number oe methods to process
   */
  int loopSize;

  /** \brief get number of evaluated element in the loop
   */
  inline int getNinTheLoop() const { return nInTheLoop; }

  /** \brief get loop size
   */
  inline int getLoopSize() const { return loopSize; }

  int rAnk; ///< processor rank

  int sIze; ///< number of processors in communicator

  const RefEntity_multiIndex
      *refinedEntitiesPtr; ///< container of mofem dof entities

  const RefElement_multiIndex
      *refinedFiniteElementsPtr; ///< container of mofem finite element entities

  const Problem *problemPtr; ///< raw pointer to problem

  const Field_multiIndex *fieldsPtr; ///< raw pointer to fields container

  const FieldEntity_multiIndex
      *entitiesPtr; ///< raw pointer to container of field entities

  const DofEntity_multiIndex *dofsPtr; ///< raw pointer container of dofs

  const FiniteElement_multiIndex
      *finiteElementsPtr; ///< raw pointer to container finite elements

  const EntFiniteElement_multiIndex
      *finiteElementsEntitiesPtr; ///< raw pointer to container finite elements
                                  ///< entities

  const FieldEntityEntFiniteElementAdjacencyMap_multiIndex
      *adjacenciesPtr; ///< raw pointer to container to adjacencies between dofs
                       ///< and finite elements

  inline unsigned int getFieldBitNumber(std::string field_name) const {
    auto field_it = fieldsPtr->get<FieldName_mi_tag>().find(field_name);
    if(field_it != fieldsPtr->get<FieldName_mi_tag>().end())
      return (*field_it)->getBitNumber();
    else
      return BITFEID_SIZE;
  }


  /**
   * @brief Copy data from other base method to this base method
   *
   * @param basic
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode copyBasicMethod(const BasicMethod &basic);

  /**
   * @brief Hook function for pre-processing
   */
  boost::function<MoFEMErrorCode()> preProcessHook;

  /**
   * @brief Hook function for post-processing
   */
  boost::function<MoFEMErrorCode()> postProcessHook;

  /**
   * @brief Hook function for operator
   */
  boost::function<MoFEMErrorCode()> operatorHook;

  /** \brief function is run at the beginning of loop
   *
   * It is used to zeroing matrices and vectors, calculation of shape
   * functions on reference element, preprocessing boundary conditions, etc.
   */
  virtual MoFEMErrorCode preProcess();

  /** \brief function is run for every finite element
   *
   * It is used to calculate element local matrices and assembly. It can be
   * used for post-processing.
   */
  virtual MoFEMErrorCode operator()();

  /** \brief function is run at the end of loop
   *
   * It is used to assembly matrices and vectors, calculating global variables,
   * f.e. total internal energy, ect.
   *
   * Iterating over dofs:
   * Example1 iterating over dofs in row by name of the field
   * for(_IT_GET_FEROW_BY_NAME_DOFS_FOR_LOOP_(this,"DISPLACEMENT",it)) { ... }
   *
   *
   */
  virtual MoFEMErrorCode postProcess();

  boost::movelib::unique_ptr<bool> vecAssembleSwitch;
  boost::movelib::unique_ptr<bool> matAssembleSwitch;
};

/**
 * \brief structure for User Loop Methods on finite elements
 * \ingroup mofem_loops
 *
 * It can be used to calculate stiffness matrices, residuals, load vectors etc.
 * It is low level class however in some class users looking for speed and
 * efficiency, can use it directly.
 *
 * This class is used with Interface::loop_finite_elements, where
 * user overloaded operator FEMethod::operator() is executed for each element in
 * the problem. Class have to additional methods which are overloaded by user,
 * FEMethod::preProcess() and FEMethod::postProcess() executed at beginning and
 * end of the loop over problem elements, respectively.
 *
 */
struct FEMethod : public BasicMethod {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const {
    MoFEMFunctionBeginHot;
    if (uuid == IDD_MOFEMFEMethod) {
      *iface = const_cast<FEMethod *>(this);
      MoFEMFunctionReturnHot(0);
    }

    ierr = query_interface(uuid, iface);
    CHKERRG(ierr);
    MoFEMFunctionReturnHot(0);
  }

  FEMethod();

  std::string feName; ///< Name of finite element

  boost::shared_ptr<const NumeredEntFiniteElement>
      numeredEntFiniteElementPtr; ///< Pointer to finite element database
                                  ///< structure

  inline const FEDofEntity_multiIndex &getDataDofs() const {
    return numeredEntFiniteElementPtr->getDataDofs();
  };

  inline boost::shared_ptr<FEDofEntity_multiIndex> &getDataDofsPtr() const {
    return const_cast<NumeredEntFiniteElement *>(
               numeredEntFiniteElementPtr.get())
        ->getDataDofsPtr();
  };

  inline const FieldEntity_multiIndex_spaceType_view &getDataFieldEnts() const {
    return numeredEntFiniteElementPtr->getDataFieldEnts();
  }

  inline boost::shared_ptr<FieldEntity_multiIndex_spaceType_view> &
  getDataFieldEntsPtr() const {
    return const_cast<NumeredEntFiniteElement *>(
               numeredEntFiniteElementPtr.get())
        ->getDataFieldEntsPtr();
  }

  inline const FieldEntity_vector_view &getRowFieldEnts() const {
    return numeredEntFiniteElementPtr->getRowFieldEnts();
  };

  inline boost::shared_ptr<FieldEntity_vector_view> &
  getRowFieldEntsPtr() const {
    return const_cast<NumeredEntFiniteElement *>(
               numeredEntFiniteElementPtr.get())
        ->getRowFieldEntsPtr();
  };

  inline const FieldEntity_vector_view &getColFieldEnts() const {
    return numeredEntFiniteElementPtr->getColFieldEnts();
  };

  inline boost::shared_ptr<FieldEntity_vector_view> &
  getColFieldEntsPtr() const {
    return const_cast<NumeredEntFiniteElement *>(
               numeredEntFiniteElementPtr.get())
        ->getColFieldEntsPtr();
  };

  inline const FENumeredDofEntity_multiIndex &getRowDofs() const {
    return numeredEntFiniteElementPtr->getRowDofs();
  };

  inline boost::shared_ptr<FENumeredDofEntity_multiIndex> &
  getRowDofsPtr() const {
    return const_cast<NumeredEntFiniteElement *>(
               numeredEntFiniteElementPtr.get())
        ->getRowDofsPtr();
  };

  inline const FENumeredDofEntity_multiIndex &getColDofs() const {
    return numeredEntFiniteElementPtr->getColDofs();
  };

  inline boost::shared_ptr<FENumeredDofEntity_multiIndex> &
  getColDofsPtr() const {
    return const_cast<NumeredEntFiniteElement *>(
               numeredEntFiniteElementPtr.get())
        ->getColDofsPtr();
  };

  /// \brief Get number of DOFs on element
  MoFEMErrorCode getNumberOfNodes(int &num_nodes) const;

  inline EntityHandle getFEEntityHandle() const;

  MoFEMErrorCode getNodeData(const std::string field_name, VectorDouble &data,
                             const bool reset_dofs = true);

  template <class MULTIINDEX>
  typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,
                                          const std::string &field_name,
                                          const EntityType type) const {
    return index.lower_bound(boost::make_tuple(field_name, type));
  }

  template <class MULTIINDEX>
  typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,
                                        const std::string &field_name,
                                        const EntityType type) const {
    return index.upper_bound(boost::make_tuple(field_name, type));
  }

/** \brief loop over all dofs which are on a particular FE data, field and
 * entity type \ingroup mofem_loops
 */
#define _IT_GET_FEDATA_BY_TYPE_DOFS_FOR_LOOP_(FE, NAME, TYPE, IT)              \
  auto IT = FE->get_begin<FEDofEntityByNameAndEnt>(                            \
      FE->dataPtr->get<FEDofEntityByNameAndEnt>(), NAME,                       \
      get_id_for_min_type(TYPE));                                              \
  IT != FE->get_end<FEDofEntityByNameAndEnt>(                                  \
            FE->dataPtr->get<FEDofEntityByNameAndEnt>(), NAME,                 \
            get_id_for_max_type(TYPE));                                        \
  IT++

  template <class MULTIINDEX>
  typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,
                                          const std::string &field_name) const {
    return index.lower_bound(field_name);
  }
  template <class MULTIINDEX>
  typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,
                                        const std::string &field_name) const {
    return index.upper_bound(field_name);
  }

/** \brief loop over all dofs which are on a particular FE row and field
 * \ingroup mofem_loops
 */
#define _IT_GET_FEROW_BY_NAME_DOFS_FOR_LOOP_(FE, NAME, IT)                     \
  auto IT = FE->get_begin<FENumeredDofEntityByUId>(                            \
      FE->rowPtr->get<Unique_mi_tag>(),                                        \
      FieldEntity::getLoBitNumberUId(getFieldBitNumber(NAME)));                \
  IT != FE->get_end<FENumeredDofEntityByUId>(                                  \
            FE->rowPtr->get<Unique_mi_tag>(),                                  \
            FieldEntity::getHiBitNumberUId(getFieldBitNumber(NAME)));          \
  IT++

/** \brief loop over all dofs which are on a particular FE column and field
 * \ingroup mofem_loops
 */
#define _IT_GET_FECOL_BY_NAME_DOFS_FOR_LOOP_(FE, NAME, IT)                     \
  auto IT = FE->get_begin<FENumeredDofEntityByUId>(                            \
      FE->colPtr->get<Unique_mi_tag>(),                                        \
      FieldEntity::getLoBitNumberUId(getFieldBitNumber(NAME)));                \
  IT != FE->get_end<FENumeredDofEntityByUId>(                                  \
            FE->colPtr->get<Unique_mi_tag>(),                                  \
            FieldEntity::getHiBitNumberUId(getFieldBitNumber(NAME)));          \
  IT++

/** \brief loop over all dofs which are on a particular FE data and field
 * \ingroup mofem_loops
 */
#define _IT_GET_FEDATA_BY_NAME_DOFS_FOR_LOOP_(FE, NAME, IT)                    \
  auto IT = FE->get_begin<FEDofEntityByUId>(                                   \
      FE->dataPtr->get<Unique_mi_tag>(),                                       \
      FieldEntity::getLoBitNumberUId(getFieldBitNumber(NAME)));                \
  IT != FE->get_end<FEDofEntityByUId>(                                         \
            FE->dataPtr->get<Unique_mi_tag>(),                                 \
            FieldEntity::getHiBitNumberUId(getFieldBitNumber(NAME)));          \
  IT++
};

inline EntityHandle FEMethod::getFEEntityHandle() const {
  return numeredEntFiniteElementPtr->getEnt();
}

/**
 * \brief Data structure to exchange data between mofem and User Loop Methods on
 * entities. \ingroup mofem_loops
 *
 * It allows to exchange data between MoFEM and user functions. It stores
 * information about multi-indices.
 */
struct EntityMethod : public BasicMethod {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const {
    MoFEMFunctionBegin;
    if (uuid == IDD_MOFEMEntityMethod) {
      *iface = const_cast<EntityMethod *>(this);
      MoFEMFunctionReturnHot(0);
    }
    CHKERR query_interface(uuid, iface);
    MoFEMFunctionReturn(0);
  }

  EntityMethod();

  boost::shared_ptr<Field> fieldPtr;
  boost::shared_ptr<FieldEntity> entPtr;
};

/**
 * \brief Data structure to exchange data between mofem and User Loop Methods on
 * entities. \ingroup mofem_loops
 *
 * It allows to exchange data between MoFEM and user functions. It stores
 * information about multi-indices.
 */
struct DofMethod : public BasicMethod {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const {
    MoFEMFunctionBeginHot;
    if (uuid == IDD_MOFEMDofMethod) {
      *iface = const_cast<DofMethod *>(this);
      MoFEMFunctionReturnHot(0);
    }

    CHKERR query_interface(uuid, iface);
    MoFEMFunctionReturnHot(0);
  }

  DofMethod();

  boost::shared_ptr<Field> fieldPtr;
  boost::shared_ptr<DofEntity> dofPtr;
  boost::shared_ptr<NumeredDofEntity> dofNumeredPtr;
};

/// \deprecated name changed use DofMethod insead EntMethod
DEPRECATED typedef DofMethod EntMethod;

} // namespace MoFEM

#endif // __LOOPMETHODS_HPP__

/**
 * \defgroup mofem_loops Loops
 * \ingroup mofem
 */
