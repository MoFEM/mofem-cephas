/** \file LoopMethods.hpp
 * \brief MoFEM interface
 *
 * Data structures for making loops over finite elements and entities in the
 * problem or MoFEM database.
 *
 */

#ifndef __LOOPMETHODS_HPP__
#define __LOOPMETHODS_HPP__

namespace MoFEM {
struct PetscData : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  PetscData();

  virtual ~PetscData() = default;

  enum DataContext {
    CTX_SET_NONE = 0,
    CTX_SET_F = 1 << 0,
    CTX_SET_A = 1 << 1,
    CTX_SET_B = 1 << 2,
    CTX_SET_X = 1 << 3,
    CTX_SET_DX = 1 << 4,
    CTX_SET_X_T = 1 << 5,
    CTX_SET_X_TT = 1 << 6,
    CTX_SET_TIME = 1 << 7
  };

  using Switches = std::bitset<8>;

  static constexpr Switches CtxSetNone = PetscData::Switches(CTX_SET_NONE);
  static constexpr Switches CtxSetF = PetscData::Switches(CTX_SET_F);
  static constexpr Switches CtxSetA = PetscData::Switches(CTX_SET_A);
  static constexpr Switches CtxSetB = PetscData::Switches(CTX_SET_B);
  static constexpr Switches CtxSetX = PetscData::Switches(CTX_SET_X);
  static constexpr Switches CtxSetDX = PetscData::Switches(CTX_SET_DX);
  static constexpr Switches CtxSetX_T = PetscData::Switches(CTX_SET_X_T);
  static constexpr Switches CtxSetX_TT = PetscData::Switches(CTX_SET_X_TT);
  static constexpr Switches CtxSetTime = PetscData::Switches(CTX_SET_TIME);

  Switches data_ctx;

  MoFEMErrorCode copyPetscData(const PetscData &petsc_data);

  Vec f;
  Mat A;
  Mat B;
  Vec x;
  Vec dx;
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

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
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

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  enum SNESContext { CTX_SNESSETFUNCTION, CTX_SNESSETJACOBIAN, CTX_SNESNONE };

  SnesMethod();

  virtual ~SnesMethod() = default;

  /**
   * \brief Copy snes data
   */
  MoFEMErrorCode copySnes(const SnesMethod &snes);

  SNESContext snes_ctx;

  SNES snes;    ///< snes solver
  Vec &snes_x;  ///< state vector
  Vec &snes_dx; ///< solution update
  Vec &snes_f;  ///< residual
  Mat &snes_A;  ///< jacobian matrix
  Mat &snes_B;  ///< preconditioner of jacobian matrix
};

/**
 * \brief data structure for TS (time stepping) context
 * \ingroup mofem_loops
 *
 * Structure stores context data which are set in functions run by PETSc Time
 * Stepping functions.
 */
struct TSMethod : virtual protected PetscData {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
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

  PetscInt ts_step; ///< time step number
  PetscReal ts_a;   ///< shift for U_t (see PETSc Time Solver)
  PetscReal ts_aa;  ///< shift for U_tt shift for U_tt
  PetscReal ts_t;   ///< time
  PetscReal ts_dt;  ///< time step size

  Vec &ts_u;    ///< state vector
  Vec &ts_u_t;  ///< time derivative of state vector
  Vec &ts_u_tt; ///< second time derivative of state vector
  Vec &ts_F;    ///< residual vector

  Mat &ts_A; ///< Jacobian of G(U) = F(t,U,W+v*U,W'+a*U), equivalent to dF/dU +
             ///< v*dF/dU_t + a*dF/dU_tt
  Mat &ts_B; ///< Preconditioner for ts_A
};

/**
 * \brief Data structure to exchange data between mofem and User Loop Methods.
 * \ingroup mofem_loops
 *
 * It allows to exchange data between MoFEM and user functions. It stores
 * information about multi-indices.
 *
 */
struct BasicMethod : public KspMethod, SnesMethod, TSMethod {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const {
    MoFEMFunctionBeginHot;
    *iface = const_cast<BasicMethod *>(this);
    MoFEMFunctionReturnHot(0);
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

  /**
   * @brief Llo and hi processor rank of iterated entities
   * 
   */
  std::pair<int, int> loHiFERank;

  /**
   * @brief Get lo and hi processor rank of iterated entities
   * 
   * @return raturn std::pair<int, int> loHiFERank
   */
  inline auto getLoHiFERank() const { return loHiFERank; }

  /**
   * @brief Get upper rank in loop for iterating elements
   * 
   * @return loHiFERank.first
   */
  inline auto getLoFERank() const { return loHiFERank.first; }

  /**
   * @brief Get upper rank in loop for iterating elements
   * 
   * @return loHiFERank.first
   */
  inline auto getHiFERank() const { return loHiFERank.second; }

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
    if (fieldsPtr) {
      auto field_it = fieldsPtr->get<FieldName_mi_tag>().find(field_name);
      if (field_it != fieldsPtr->get<FieldName_mi_tag>().end())
        return (*field_it)->getBitNumber();
      else
        return BITFEID_SIZE;
    } else {
      THROW_MESSAGE("Pointer to fields multi-index is not set");
      return BITFEID_SIZE;
    }
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

  /**
   * @brief Get the cache weak ptr object
   *
   * \note This store problem information on entities about DOFs. Each problem
   * store different information. If you iterate over finite elements in
   * preprocessor of TS solve element, us TS cache in the loop. Otherwise you
   * will create undetermined behaviour or segmentation error. This is necessary
   * compromise over bug resilience for memory saving and performance.
   *
   * @return boost::weak_ptr<CacheTuple>
   */
  inline boost::weak_ptr<CacheTuple> getCacheWeakPtr() const {
    return cacheWeakPtr;
  }

  boost::movelib::unique_ptr<bool> vecAssembleSwitch;
  boost::movelib::unique_ptr<bool> matAssembleSwitch;

  boost::weak_ptr<CacheTuple> cacheWeakPtr; // cache pointer entity data
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

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const {
    MoFEMFunctionBeginHot;
    *iface = const_cast<FEMethod *>(this);
    MoFEMFunctionReturnHot(0);
  }

  FEMethod() = default;

  std::string feName; ///< Name of finite element

  boost::shared_ptr<const NumeredEntFiniteElement>
      numeredEntFiniteElementPtr; ///< Pointer to finite element database
                                  ///< structure

  /**
   * @brief get finite element name
   *
   * @return std::string
   */
  inline auto getFEName() const { return feName; }

  /**
   * @brief Tet if element to skip element
   *
   * If is set and return false  elemnent us skiped in
   * MoFEM::Core::loop_finite_elements
   * 
   * \note That functionality is used to run elements on particular bit levels
   *
   */
  boost::function<bool(FEMethod *fe_method_ptr)> exeTestHook;

  inline auto getDataDofsPtr() const {
    return numeredEntFiniteElementPtr->getDataDofsPtr();
  };

  inline auto getDataVectorDofsPtr() const {
    return numeredEntFiniteElementPtr->getDataVectorDofsPtr();
  };

  inline const FieldEntity_vector_view &getDataFieldEnts() const {
    return numeredEntFiniteElementPtr->getDataFieldEnts();
  }

  inline boost::shared_ptr<FieldEntity_vector_view> &
  getDataFieldEntsPtr() const {
    return const_cast<NumeredEntFiniteElement *>(
               numeredEntFiniteElementPtr.get())
        ->getDataFieldEntsPtr();
  }

  inline auto &getRowFieldEnts() const {
    return numeredEntFiniteElementPtr->getRowFieldEnts();
  };

  inline auto &getRowFieldEntsPtr() const {
    return numeredEntFiniteElementPtr->getRowFieldEntsPtr();
  };

  inline auto &getColFieldEnts() const {
    return numeredEntFiniteElementPtr->getColFieldEnts();
  };

  inline auto &getColFieldEntsPtr() const {
    return numeredEntFiniteElementPtr->getColFieldEntsPtr();
  };

  inline auto getRowDofsPtr() const {
    return numeredEntFiniteElementPtr->getRowDofsPtr();
  };

  inline auto getColDofsPtr() const {
    return numeredEntFiniteElementPtr->getColDofsPtr();
  };

  inline auto getNumberOfNodes() const;

  inline EntityHandle getFEEntityHandle() const;

  MoFEMErrorCode getNodeData(const std::string field_name, VectorDouble &data,
                             const bool reset_dofs = true);

};

inline auto FEMethod::getNumberOfNodes() const {
  return moab::CN::VerticesPerEntity(numeredEntFiniteElementPtr->getEntType());
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

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const {
    MoFEMFunctionBegin;
    *iface = const_cast<EntityMethod *>(this);
    MoFEMFunctionReturn(0);
  }

  EntityMethod() = default;

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

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const {
    MoFEMFunctionBeginHot;
    *iface = const_cast<DofMethod *>(this);
    MoFEMFunctionReturnHot(0);
  }

  DofMethod() = default;

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
