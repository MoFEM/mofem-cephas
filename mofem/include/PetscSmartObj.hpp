/** \file PetscSmartObj.hpp
 * \brief Petsc smart obj declarations
 */

#ifndef __PETSCSPARTOBJ_HPP__
#define __PETSCSPARTOBJ_HPP__

namespace MoFEM {
template <typename T> inline PetscObject getPetscObject(T obj) {
  return reinterpret_cast<PetscObject>(obj);
}
} // namespace MoFEM

/**
 * @brief It is used by intrusive_ptr to bump reference
 *
 * \note It should not be used directly, it is internally called by
 * intrusive_ptr
 *
 * @tparam OBJ
 * @param obj
 */
template <typename OBJ> void intrusive_ptr_add_ref(OBJ obj) {
  if (obj) {
    PetscErrorCode ierr = PetscObjectReference(MoFEM::getPetscObject(obj));
    CHKERRABORT(PetscObjectComm(MoFEM::getPetscObject(obj)), ierr);
  }
}

/**
 * @brief It is used by intrusive_ptr to dereference and destroy petsc object
 *
 * \note It should not be used directly, it is internally called by
 * intrusive_ptr
 *
 * @tparam OBJ
 * @param obj
 */
template <typename OBJ> void intrusive_ptr_release(OBJ obj) {
  if (obj) {
    int cnt = 0;
    PetscErrorCode ierr =
        PetscObjectGetReference(MoFEM::getPetscObject(obj), &cnt);
    if (!ierr) {
      if (cnt) {
        if (cnt > 1) {
          ierr = PetscObjectDereference(MoFEM::getPetscObject(obj));
        } else {
          ierr = PetscObjectDestroy(reinterpret_cast<PetscObject *>(&obj));
        }
      }
      auto comm = PetscObjectComm(MoFEM::getPetscObject(obj));
      CHKERRABORT(comm, ierr);
    }
  }
}

template <> void intrusive_ptr_release<Vec>(Vec obj);
template <> void intrusive_ptr_release<Mat>(Mat obj);
template <> void intrusive_ptr_release<DM>(DM obj);
template <> void intrusive_ptr_release<IS>(IS obj);
template <> void intrusive_ptr_release<AO>(AO obj);
template <> void intrusive_ptr_release<KSP>(KSP obj);
template <> void intrusive_ptr_release<SNES>(SNES obj);
template <> void intrusive_ptr_release<TS>(TS obj);

namespace MoFEM {
/**
 * @brief intrusive_ptr for managing petsc objects
 *
 * It manages destruction, referencing and dereferencing petsc objects. It is
 * similar how smart_ptr pointers works, but applied for petsc objects like Vec,
 * DM, Mat, etc.
 *
 * \code
 * SmartPetscObj<Vec> smart_vec = createSmartGhostVector(...);
 * \endcode
 *
 * @tparam OBJ
 */
template <typename OBJ>
struct SmartPetscObj
    : public boost::intrusive_ptr<typename std::remove_pointer<OBJ>::type> {

  using Derived = boost::intrusive_ptr<typename std::remove_pointer<OBJ>::type>;

  using Derived::Derived;

  SmartPetscObj(std::nullptr_t ptr) : SmartPetscObj() {}

  /**
   * @brief Construct a new Smart Petsc Obj object
   *
   * \note If add_red is set to true, you have to destroy OBJ.
   *
   * @param o
   * @param add_ref // if false ownership of OBJ is taken by SmartPetscObj
   */
  explicit SmartPetscObj(OBJ o, bool add_ref = false)
      : boost::intrusive_ptr<typename std::remove_pointer<OBJ>::type>(o,
                                                                      add_ref) {
  }

  operator OBJ() { return this->get(); }
  explicit operator PetscObject() {
    return reinterpret_cast<PetscObject>(this->get());
  }

  int use_count() const {
    if (this->get()) {
      int cnt;
      ierr = PetscObjectGetReference(getPetscObject(this->get()), &cnt);
      CHKERRABORT(PetscObjectComm(getPetscObject(this->get())), ierr);
      return cnt;
    } else
      return 0;
  }
};

/**
 * @brief Creates smart DM object
 *
 * DM object can be used as any other object, but is destroyed as smart pointer
 * when no longer used.
 *
 * \code
 * CHKERR DMRegister_MoFEM("MOFEM")
 * {
 *    auto dm = createDM(PETSC_COMM_WORLD, "MOFEM");
 *
 *    // ...
 *
 *    // dm is automatically destroyed when program goes out of the scope
 * }
 *
 *
 *
 * \endcode
 *
 */
inline auto createDM(MPI_Comm comm, const std::string dm_type_name) {
  DM dm;
  CHK_THROW_MESSAGE(DMCreate(comm, &dm), "Failed to create DM");
  CHK_THROW_MESSAGE(DMSetType(dm, dm_type_name.c_str()), "Failed set DM type");
  return SmartPetscObj<DM>(dm);
};

/** @deprecated use createDM */
DEPRECATED inline auto createSmartDM(MPI_Comm comm,
                                     const std::string dm_type_name) {
  return createDM(comm, dm_type_name);
}

/**
 * @brief Get the Comm From Petsc Object object
 *
 * @param obj
 * @return MPI_Comm
 */
inline MPI_Comm getCommFromPetscObject(PetscObject obj) {
  MPI_Comm comm;
  CHK_THROW_MESSAGE(PetscObjectGetComm(obj, &comm),
                    "Failed to get comm from PETSc object");
  return comm;
};

/**
 * @brief Create smart ghost vector
 *
 * For details abut arguments see here:
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecCreateGhost.html>VecCreateGhost</a>.
 *
 * \code
 * auto vec = createGhostVector(...);
 * \endcode
 *
 */
inline auto createGhostVector(MPI_Comm comm, PetscInt n, PetscInt N,
                              PetscInt nghost, const PetscInt ghosts[]) {
  Vec vv;
  CHK_THROW_MESSAGE(VecCreateGhost(comm, n, N, nghost, ghosts, &vv),
                    "Failed to create ghosted Vec");
  return SmartPetscObj<Vec>(vv);
};

/** @deprecated use createGhostVector */
DEPRECATED inline auto createSmartGhostVector(MPI_Comm comm, PetscInt n,
                                              PetscInt N, PetscInt nghost,
                                              const PetscInt ghosts[]) {
  return createGhostVector(comm, n, N, nghost, ghosts);
}

/**
 * @brief Create MPI Vector
 *
 * For details abut arguments see here:
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecCreateMPI.html>VecCreateMPI</a>.
 *
 */
inline auto createVectorMPI(MPI_Comm comm, PetscInt n, PetscInt N) {
  Vec vv;
  CHK_THROW_MESSAGE(VecCreateMPI(comm, n, N, &vv), "Failed to create Vec");
  return SmartPetscObj<Vec>(vv);
};

/** @deprecated use createVectorMPI */
DEPRECATED inline auto createSmartVectorMPI(MPI_Comm comm, PetscInt n,
                                            PetscInt N) {
  return createVectorMPI(comm, n, N);
}

/**
 * @brief Create duplicate vector of smart vector
 *
 * For details abut arguments see here:
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecDuplicate.html>VecDuplicate</a>.
 */
inline SmartPetscObj<Vec> vectorDuplicate(Vec vec) {
  Vec duplicate;
  CHK_THROW_MESSAGE(VecDuplicate(vec, &duplicate), "Failed to duplicate Vec");
  return SmartPetscObj<Vec>(duplicate);
};

/**
 * @deprecated use vectorDuplicate
 */
DEPRECATED inline SmartPetscObj<Vec> smartVectorDuplicate(Vec vec) {
  return vectorDuplicate(vec);
}

inline SmartPetscObj<Mat> matDuplicate(Mat mat, MatDuplicateOption op) {
  Mat duplicate;
  CHK_THROW_MESSAGE(MatDuplicate(mat, op, &duplicate),
                    "Failed to duplicate Mat");
  return SmartPetscObj<Mat>(duplicate);
};

/**
 * @deprecated use matDuplicate
 */
DEPRECATED inline SmartPetscObj<Mat> smartMatDuplicate(Mat mat,
                                                       MatDuplicateOption op) {
  return matDuplicate(mat, op);
}

inline auto createTS(MPI_Comm comm) {
  TS ts;
  CHK_THROW_MESSAGE(TSCreate(comm, &ts), "Failed to create TS");
  return SmartPetscObj<TS>(ts);
};

inline auto createSNES(MPI_Comm comm) {
  SNES snes;
  CHK_THROW_MESSAGE(SNESCreate(comm, &snes), "Failed to create SNES");
  return SmartPetscObj<SNES>(snes);
};

inline auto createKSP(MPI_Comm comm) {
  KSP ksp;
  CHK_THROW_MESSAGE(KSPCreate(comm, &ksp), "Failed to create KSP");
  return SmartPetscObj<KSP>(ksp);
};

inline auto createPC(MPI_Comm comm) {
  PC pc;
  CHK_THROW_MESSAGE(PCCreate(comm, &pc), "Failed to create PC");
  return SmartPetscObj<PC>(pc);
};

/**
 * @brief Creates a data structure for an index set containing a list of
 * integers.
 *
 * <a
 * href=https://petsc.org/release/docs/manualpages/IS/ISCreateGeneral/>AOCreateMappingIS</a>.
 *
 * @param comm the MPI communicator
 * @param n the length of the index set
 * @param idx  the list of integers
 * @param mode PETSC_COPY_VALUES, PETSC_OWN_POINTER, or PETSC_USE_POINTER; see
 * PetscCopyMode for meaning of this flag.
 * @return SmartPetscObj<IS>(is)
 */
inline auto createISGeneral(MPI_Comm comm, PetscInt n, const PetscInt idx[],
                            PetscCopyMode mode) {
  IS is;
  CHK_THROW_MESSAGE(ISCreateGeneral(comm, n, idx, mode, &is), "Create IS");
  return SmartPetscObj<IS>(is);
}

/**
 * @brief IS All gather
 * 
 * @param is 
 * @return auto 
 */
inline auto  isAllGather(IS is) {
  IS isout;
  CHK_THROW_MESSAGE(ISAllGather(is, &isout), "Failed to create ISAllGather");
  return SmartPetscObj<IS>(isout);
}

/**
 * @brief Creates an application mapping using two index sets.
 *
 * <a
 * href=https://petsc.org/release/docs/manualpages/AO/AOCreateMappingIS/>AOCreateMappingIS</a>.
 *
 * @param isapp  index set that defines an ordering
 * @param ispetsc  index set that defines another ordering, maybe NULL for
 * identity
 * @param aoout the new application ordering
 * @return SmartPetscObj<AO>(ao)
 */
inline auto createAOMappingIS(IS isapp, IS ispetsc) {
  AO ao;
  CHK_THROW_MESSAGE(AOCreateMappingIS(isapp, ispetsc, &ao),
                    "Failed to create AO");
  return SmartPetscObj<AO>(ao);
};

/**
 * @brief Creates an application mapping using two integer arrays.
 *
 * <a
 * href=https://petsc.org/release/docs/manualpages/AO/AOCreateMapping/>AOCreateMappingIS</a>.
 *
 * @param comm MPI communicator that is to share the AO
 * @param napp size of integer arrays
 * @param myapp  integer array that defines an ordering
 * @param mypetsc  integer array that defines another ordering (may be NULL to
 * indicate the identity ordering)
 * @return SmartPetscObj<AO>(ao);
 */
inline auto createAOMapping(MPI_Comm comm, PetscInt napp,
                            const PetscInt myapp[], const PetscInt mypetsc[]) {
  AO ao;
  CHK_THROW_MESSAGE(AOCreateMapping(comm, napp, myapp, mypetsc, &ao),
                    "create ao");
  return SmartPetscObj<AO>(ao);
}

/**
 * @brief Create a Vec Scatter object
 *
 * <a
 * href=https://petsc.org/release/manualpages/PetscSF/VecScatterCreate/>VecScatterCreate</a>.
 *
 * @param x a vector that defines the shape (parallel data layout of the vector)
 * of vectors from which we scatter
 * @param ix the indices of xin to scatter (if NULL scatters all values)
 * @param y a vector that defines the shape (parallel data layout of the vector)
 * of vectors to which we scatter
 * @param iy the indices of yin to hold results (if NULL fills entire vector yin
 * in order)
 * @return
 */
inline auto createVecScatter(Vec x, IS ix, Vec y, IS iy) {
  VecScatter s;
  CHK_THROW_MESSAGE(VecScatterCreate(x, ix, y, iy, &s), "create scatter");
  return SmartPetscObj<VecScatter>(s);
}

/**
 * @brief Get ISDifference
 *
 * <a
 * href=https://petsc.org/release/docs/manualpages/IS/ISDifference/>ISDifference</a>.
 *
 * @param is1 first index, to have items removed from it
 * @param is2 index values to be removed
 * @return is1 - is2
 */
inline auto isDifference(IS is1, IS is2) {
  IS is_raw;
  CHK_THROW_MESSAGE(ISDifference(is1, is2, &is_raw), "create difference");
  return SmartPetscObj<IS>(is_raw);
}

inline auto createISLocalToGlobalMapping(IS is) {
  ISLocalToGlobalMapping map_raw;
  CHK_THROW_MESSAGE(ISLocalToGlobalMappingCreateIS(is, &map_raw),
                    "create local to global mapping");
  return SmartPetscObj<ISLocalToGlobalMapping>(map_raw);
}

} // namespace MoFEM

#endif