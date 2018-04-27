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

/**
 * \brief data structure for ksp (linear solver) context
 * \ingroup mofem_loops
 *
 * Struture stores context data which are set in functions run by PETSc SNES
 * functions.
 *
 */
struct KspMethod : virtual public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const {
    MoFEMFunctionBeginHot;
    if (uuid == IDD_MOFEMKspMethod) {
      *iface = const_cast<KspMethod *>(this);
      MoFEMFunctionReturnHot(0);
    }
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  }

  /**
   * \brief pass information about context of KSP/DM for with finite element is
   * computed
   */
  enum KSPContext { CTX_SETFUNCTION, CTX_OPERATORS, CTX_KSPNONE };

  KspMethod()
      : ksp_ctx(CTX_KSPNONE), ksp(PETSC_NULL), ksp_f(PETSC_NULL),
        ksp_A(PETSC_NULL), ksp_B(PETSC_NULL) {}

  virtual ~KspMethod() {}

  /**
   * \brief set operator type
   * @param  ctx Context, CTX_SETFUNCTION, CTX_OPERATORS, CTX_KSPNONE
   * @return     error code
   */
  MoFEMErrorCode setKspCtx(const KSPContext &ctx);

  /**
   * \brief set solver
   * @param  ksp solver
   * @return     error code
   */
  MoFEMErrorCode setKsp(KSP ksp);

  /**
   * \brief copy data form another method
   * @param  ksp ksp method
   * @return     error code
   */
  MoFEMErrorCode copyKsp(const KspMethod &ksp);

  KSPContext ksp_ctx; ///< Context
  KSP ksp;            ///< KSP solver
  Vec ksp_f;          ///< the right hand side vector
  Mat ksp_A;          ///< matrix
  Mat ksp_B;          ///< preconditioner matrix
};

/**
 * \brief data structure for snes (nonlinear solver) context
 * \ingroup mofem_loops
 *
 * Structure stores context data which are set in functions run by PETSc SNES
 * functions.
 *
 */
struct SnesMethod : virtual public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const {
    if (uuid == IDD_MOFEMSnesMethod) {
      *iface = const_cast<SnesMethod *>(this);
      MoFEMFunctionReturnHot(0);
    }
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  }

  enum SNESContext { CTX_SNESSETFUNCTION, CTX_SNESSETJACOBIAN, CTX_SNESNONE };

  SNESContext snes_ctx;
  SnesMethod()
      : snes_ctx(CTX_SNESNONE), snes_x(PETSC_NULL), snes_f(PETSC_NULL),
        snes_A(PETSC_NULL), snes_B(PETSC_NULL) {}

  virtual ~SnesMethod() {}

  /**
   * \brief Set SNES context
   */
  MoFEMErrorCode setSnesCtx(const SNESContext &ctx);

  /**
   * \brief Set SNES instance
   */
  MoFEMErrorCode setSnes(SNES snes);

  /**
   * \brief Copy snes data
   */
  MoFEMErrorCode copySnes(const SnesMethod &snes);

  SNES snes;
  Vec snes_x, snes_f;
  Mat snes_A, snes_B;

};

/**
 * \brief data structure for TS (time stepping) context
 * \ingroup mofem_loops
 *
 * Structure stores context data which are set in functions run by PETSc Time
 * Stepping functions.
 */
struct TSMethod : virtual public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const {
    if (uuid == IDD_MOFEMTsMethod) {
      *iface = const_cast<TSMethod *>(this);
      MoFEMFunctionReturnHot(0);
    }
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  }

  enum TSContext {
    CTX_TSSETRHSFUNCTION,
    CTX_TSSETRHSJACOBIAN,
    CTX_TSSETIFUNCTION,
    CTX_TSSETIJACOBIAN,
    CTX_TSTSMONITORSET,
    CTX_TSNONE
  };

  TSContext ts_ctx;
  TSMethod()
      : ts_ctx(CTX_TSNONE), ts_u(PETSC_NULL), ts_u_t(PETSC_NULL),
        ts_F(PETSC_NULL), ts_A(PETSC_NULL), ts_B(PETSC_NULL), ts_step(-1),
        ts_a(0), ts_t(0) {}

  virtual ~TSMethod() {}

  /// \brief Set Ts context
  MoFEMErrorCode setTsCtx(const TSContext &ctx);

  /// \brief Copy TS solver data
  MoFEMErrorCode copyTs(const TSMethod &ts);

  /// \brief Set TS solver
  MoFEMErrorCode setTs(TS _ts);

  TS ts;
  Vec ts_u, ts_u_t, ts_F;
  Mat ts_A, ts_B;

  PetscInt ts_step;
  PetscReal ts_a, ts_t;

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

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const {
    if (uuid == IDD_MOFEMBasicMethod) {
      *iface = const_cast<BasicMethod *>(this);
      MoFEMFunctionReturnHot(0);
    }
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  }

  BasicMethod();
  virtual ~BasicMethod(){};

  int nInTheLoop;
  int loopSize;

  /** \brief get number of evaluated element in the loop
   */
  inline int getNinTheLoop() const { return nInTheLoop; }

  /** \brief get loop size
   */
  inline int getLoopSize() const { return loopSize; }

  int rAnk, sIze;
  const RefEntity_multiIndex *refinedEntitiesPtr;
  const RefElement_multiIndex *refinedFiniteElementsPtr;
  const Problem *problemPtr;
  const Field_multiIndex *fieldsPtr;
  const FieldEntity_multiIndex *entitiesPtr;
  const DofEntity_multiIndex *dofsPtr;
  const FiniteElement_multiIndex *finiteElementsPtr;
  const EntFiniteElement_multiIndex *finiteElementsEntitiesPtr;
  const FieldEntityEntFiniteElementAdjacencyMap_multiIndex *adjacenciesPtr;

  MoFEMErrorCode copyBasicMethod(const BasicMethod &basic);

  /**
   * Hook function for pre-processing
   */
  boost::function<MoFEMErrorCode()> preProcessHook;

  /**
   * Hook function for post-processing
   */
  boost::function<MoFEMErrorCode()> postProcessHook;

  /**
   * Hook function for operator
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

private:
  void iNit();
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

  std::string feName;

  boost::shared_ptr<const NumeredEntFiniteElement> numeredEntFiniteElementPtr;
  boost::shared_ptr<const FENumeredDofEntity_multiIndex> rowPtr;
  boost::shared_ptr<const FENumeredDofEntity_multiIndex> colPtr;
  boost::shared_ptr<const FEDofEntity_multiIndex> dataPtr;

/** \brief loop over all dofs which are on a particular FE row
 * \ingroup mofem_loops
 */
#define _IT_GET_FEROW_DOFS_FOR_LOOP_(FE, IT)                                   \
  auto IT = FE->rowPtr->begin();                                               \
  IT != FE->rowPtr->end();                                                     \
  IT++

/** \brief loop over all dofs which are on a particular FE column
 * \ingroup mofem_loops
 */
#define _IT_GET_FECOL_DOFS_FOR_LOOP_(FE, IT)                                   \
  auto IT = FE->colPtr->begin();                                               \
  IT != FE->colPtr->end();                                                     \
  IT++

/** \brief loop over all dofs which are on a particular FE data
 * \ingroup mofem_loops
 */
#define _IT_GET_FEDATA_DOFS_FOR_LOOP_(FE, IT)                                  \
  auto IT = FE->dataPtr->begin();                                              \
  IT != FE->dataPtr->end();                                                    \
  IT++

  template <class MULTIINDEX>
  typename MULTIINDEX::iterator
  get_begin(const MULTIINDEX &index, const std::string &field_name,
            const EntityType type, const int side_number) const {
    return index.lower_bound(boost::make_tuple(field_name, type, side_number));
  }
  template <class MULTIINDEX>
  typename MULTIINDEX::iterator
  get_end(const MULTIINDEX &index, const std::string &field_name,
          const EntityType type, const int side_number) const {
    return index.upper_bound(boost::make_tuple(field_name, type, side_number));
  }

    /** \brief loop over all dofs which are on a particular FE row, field,
     * entity type and canonical side number \ingroup mofem_loops
     *
     * \param FE finite elements
     * \param Name field name
     * \param Type moab entity type (MBVERTEX, MBEDGE etc)
     * \param Side side canonical number
     * \param IT the interator in use
     */
#define _IT_GET_FEROW_BY_SIDE_DOFS_FOR_LOOP_(FE, NAME, TYPE, SIDE, IT)         \
  FENumeredDofEntity_multiIndex::index<Composite_mi_tag>::type::iterator IT =  \
      FE->get_begin<                                                           \
          FENumeredDofEntity_multiIndex::index<Composite_mi_tag>::type>(       \
          FE->rowPtr->get<Composite_mi_tag>(), NAME, TYPE, SIDE);              \
  IT != FE->get_end<                                                           \
            FENumeredDofEntity_multiIndex::index<Composite_mi_tag>::type>(     \
            FE->rowPtr->get<Composite_mi_tag>(), NAME, TYPE, SIDE);            \
  IT++

/** \brief loop over all dofs which are on a particular FE column, field, entity
 * type and canonical side number \ingroup mofem_loops
 */
#define _IT_GET_FECOL_BY_SIDE_DOFS_FOR_LOOP_(FE, NAME, TYPE, SIDE, IT)         \
  FENumeredDofEntity_multiIndex::index<Composite_mi_tag>::type::iterator IT =  \
      FE->get_begin<                                                           \
          FENumeredDofEntity_multiIndex::index<Composite_mi_tag>::type>(       \
          FE->colPtr->get<Composite_mi_tag>(), NAME, TYPE, SIDE);              \
  IT != FE->get_end<                                                           \
            FENumeredDofEntity_multiIndex::index<Composite_mi_tag>::type>(     \
            FE->colPtr->get<Composite_mi_tag>(), NAME, TYPE, SIDE);            \
  IT++

/** \brief loop over all dofs which are on a particular FE data, field, entity
 * type and canonical side number \ingroup mofem_loops
 */
#define _IT_GET_FEDATA_BY_SIDE_DOFS_FOR_LOOP_(FE, NAME, TYPE, SIDE, IT)        \
  FEDofEntity_multiIndex::index<Composite_mi_tag>::type::iterator IT =         \
      FE->get_begin<FEDofEntity_multiIndex::index<Composite_mi_tag>::type>(    \
          FE->dataPtr->get<Composite_mi_tag>(), NAME, TYPE, SIDE);             \
  IT != FE->get_end<FEDofEntity_multiIndex::index<Composite_mi_tag>::type>(    \
            FE->dataPtr->get<Composite_mi_tag>(), NAME, TYPE, SIDE);           \
  IT++

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

/** \brief loop over all dofs which are on a particular FE row, field and entity
 * type \ingroup mofem_loops
 */
#define _IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(FE, NAME, TYPE, IT)               \
  auto IT = FE->get_begin<FENumeredDofEntityByNameAndType>(                    \
      FE->rowPtr->get<Composite_Name_And_Type_mi_tag>(), NAME, TYPE);          \
  IT != FE->get_end<FENumeredDofEntityByNameAndType>(                          \
            FE->rowPtr->get<Composite_Name_And_Type_mi_tag>(), NAME, TYPE);    \
  IT++

/** \brief loop over all dofs which are on a particular FE column, field and
 * entity type \ingroup mofem_loops
 */
#define _IT_GET_FECOL_BY_TYPE_DOFS_FOR_LOOP_(FE, NAME, TYPE, IT)               \
  auto IT = FE->get_begin<FENumeredDofEntityByNameAndType>(                    \
      FE->colPtr->get<Composite_Name_And_Type_mi_tag>(), NAME, TYPE);          \
  IT != FE->get_end<FENumeredDofEntityByNameAndType>(                          \
            FE->colPtr->get<Composite_Name_And_Type_mi_tag>(), NAME, TYPE);    \
  IT++

/** \brief loop over all dofs which are on a particular FE data, field and
 * entity type \ingroup mofem_loops
 */
#define _IT_GET_FEDATA_BY_TYPE_DOFS_FOR_LOOP_(FE, NAME, TYPE, IT)              \
  auto IT = FE->get_begin<FEDofEntityByNameAndType>(                           \
      FE->dataPtr->get<Composite_Name_And_Type_mi_tag>(), NAME, TYPE);         \
  IT != FE->get_end<FEDofEntityByNameAndType>(                                 \
            FE->dataPtr->get<Composite_Name_And_Type_mi_tag>(), NAME, TYPE);   \
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
  auto IT = FE->get_begin<FENumeredDofEntityByFieldName>(                      \
      FE->rowPtr->get<FieldName_mi_tag>(), NAME);                              \
  IT != FE->get_end<FENumeredDofEntityByFieldName>(                            \
            FE->rowPtr->get<FieldName_mi_tag>(), NAME);                        \
  IT++

/** \brief loop over all dofs which are on a particular FE column and field
 * \ingroup mofem_loops
 */
#define _IT_GET_FECOL_BY_NAME_DOFS_FOR_LOOP_(FE, NAME, IT)                     \
  auto IT = FE->get_begin<FENumeredDofEntityByFieldName>(                      \
      FE->colPtr->get<FieldName_mi_tag>(), NAME);                              \
  IT != FE->get_end<FENumeredDofEntityByFieldName>(                            \
            FE->colPtr->get<FieldName_mi_tag>(), NAME);                        \
  IT++

/** \brief loop over all dofs which are on a particular FE data and field
 * \ingroup mofem_loops
 */
#define _IT_GET_FEDATA_BY_NAME_DOFS_FOR_LOOP_(FE, NAME, IT)                    \
  auto IT = FE->get_begin<FEDofEntityByFieldName>(                             \
      FE->dataPtr->get<FieldName_mi_tag>(), NAME);                             \
  IT != FE->get_end<FEDofEntityByFieldName>(                                   \
            FE->dataPtr->get<FieldName_mi_tag>(), NAME);                       \
  IT++

  template <class MULTIINDEX>
  typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,
                                          const EntityHandle ent) const {
    return index.lower_bound(ent);
  }
  template <class MULTIINDEX>
  typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,
                                        const EntityHandle ent) const {
    return index.upper_bound(ent);
  }

/** \brief loop over all dofs which are on a particular FE row and given element
 * entity (handle from moab) \ingroup mofem_loops
 */
#define _IT_GET_FEROW_DOFS_BY_ENT_FOR_LOOP_(FE, ENT, IT)                       \
  auto IT = FE->get_begin<FENumeredDofEntityByEnt>(                            \
      FE->rowPtr->get<Ent_mi_tag>(), ENT);                                     \
  IT != FE->get_end<FENumeredDofEntityByEnt>(FE->rowPtr->get<Ent_mi_tag>(),    \
                                             ENT);                             \
  IT++

/** \brief loop over all dofs which are on a particular FE column and given
 * element entity (handle from moab) \ingroup mofem_loops
 */
#define _IT_GET_FECOL_DOFS_BY_ENT_FOR_LOOP_(FE, ENT, IT)                       \
  auto IT = FE->get_begin<FENumeredDofEntityByEnt>(                            \
      FE->colPtr->get<Ent_mi_tag>(), ENT);                                     \
  IT != FE->get_end<FENumeredDofEntityByEnt>(FE->colPtr->get<Ent_mi_tag>(),    \
                                             ENT);                             \
  IT++

/** \brief loop over all dofs which are on a particular FE data and given
 * element entity (handle from moab) \ingroup mofem_loops
 */
#define _IT_GET_FEDATA_DOFS_BY_ENT_FOR_LOOP_(FE, ENT, IT)                      \
  auto IT = FE->get_begin<FEDofEntity_multiIndex::index<Ent_mi_tag>::type>(    \
      FE->dataPtr->get<Ent_mi_tag>(), ENT);                                    \
  IT != FE->get_end<FEDofEntity_multiIndex::index<Ent_mi_tag>::type>(          \
            FE->dataPtr->get<Ent_mi_tag>(), ENT);                              \
  IT++

  template <class MULTIINDEX>
  typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,
                                          const std::string &field_name,
                                          const EntityHandle ent) const {
    return index.lower_bound(boost::make_tuple(field_name, ent));
  }
  template <class MULTIINDEX>
  typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,
                                        const std::string &field_name,
                                        const EntityHandle ent) const {
    return index.upper_bound(boost::make_tuple(field_name, ent));
  }

/** \brief loop over all dofs which are on a particular FE row, field and given
 * element entity (handle from moab) \ingroup mofem_loops
 */
#define _IT_GET_FEROW_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE, NAME, ENT, IT)        \
  auto IT = FE->get_begin<FENumeredDofEntityByNameAndEnt>(                     \
      FE->rowPtr->get<Composite_Name_And_Ent_mi_tag>(), NAME, ENT);            \
  IT != FE->get_end<FENumeredDofEntityByNameAndEnt>(                           \
            FE->rowPtr->get<Composite_Name_And_Ent_mi_tag>(), NAME, ENT);      \
  IT++

/** \brief loop over all dofs which are on a particular FE column, field and
 * given element entity (handle from moab) \ingroup mofem_loops
 */
#define _IT_GET_FECOL_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE, NAME, ENT, IT)        \
  auto IT = FE->get_begin<FENumeredDofEntityByNameAndEnt>(                     \
      FE->colPtr->get<Composite_Name_And_Ent_mi_tag>(), NAME, ENT);            \
  IT != FE->get_end<FENumeredDofEntityByNameAndEnt>(                           \
            FE->colPtr->get<Composite_Name_And_Ent_mi_tag>(), NAME, ENT);      \
  IT++

/** \brief loop over all dofs which are on a particular FE data, field and given
 * element entity (handle from moab) \ingroup mofem_loops
 */
#define _IT_GET_FEDATA_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE, NAME, ENT, IT)       \
  auto IT = FE->get_begin<FEDofEntityByNameAndEnt>(                            \
      FE->dataPtr->get<Composite_Name_And_Ent_mi_tag>(), NAME, ENT);           \
  IT != FE->get_end<FEDofEntityByNameAndEnt>(                                  \
            FE->dataPtr->get<Composite_Name_And_Ent_mi_tag>(), NAME, ENT);     \
  IT++
};

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
    if (uuid == IDD_MOFEMDofMethod) {
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

/***************************************************************************/ /**
  * \defgroup mofem_loops Loops
  * \ingroup mofem
 ******************************************************************************/
