/** \file PetscSmartObj.hpp
 * \brief Petsc smart obj declarations
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
  PetscErrorCode ierr = PetscObjectReference(MoFEM::getPetscObject(obj));
  CHKERRABORT(PetscObjectComm(MoFEM::getPetscObject(obj)), ierr);
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
 * when no loneger used.
 *
 * \code
 * CHKERR DMRegister_MoFEM("MOFEM")
 * {
 *    auto dm = createSmartDM(PETSC_COMM_WORLD, "MOFEM");
 *
 *    // ...
 *
 *    // dm is autmatically destroyed when program goes out of the scope
 * }
 *
 *
 *
 * \endcode
 *
 */
inline auto createSmartDM(MPI_Comm comm, const std::string dm_type_name) {
  DM dm;
  CHK_THROW_MESSAGE(DMCreate(comm, &dm), "Failed to create DM");
  CHK_THROW_MESSAGE(DMSetType(dm, dm_type_name.c_str()), "Failed set DM type");
  return SmartPetscObj<DM>(dm);
};

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
 * auto vec = createSmartGhostVector(...);
 * \endcode
 *
 */
inline auto createSmartGhostVector(MPI_Comm comm, PetscInt n, PetscInt N,
                                 PetscInt nghost, const PetscInt ghosts[]) {
  Vec vv;
  CHK_THROW_MESSAGE(VecCreateGhost(comm, n, N, nghost, ghosts, &vv),
                    "Failed to create ghosted Vec");
  return SmartPetscObj<Vec>(vv);
};

/**
 * @brief Create MPI Vector
 *
 * For details abut arguments see here:
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecCreateMPI.html>VecCreateMPI</a>.
 *
 */
inline auto createSmartVectorMPI(MPI_Comm comm, PetscInt n, PetscInt N) {
  Vec vv;
  CHK_THROW_MESSAGE(VecCreateMPI(comm, n, N, &vv), "Failed to create Vec");
  return SmartPetscObj<Vec>(vv);
};

/**
 * @brief Create duplicate vector of smart vector
 *
 * For details abut arguments see here:
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecDuplicate.html>VecDuplicate</a>.
 */
inline SmartPetscObj<Vec> smartVectorDuplicate(SmartPetscObj<Vec> &vec) {
  if (vec.use_count()) {
    Vec duplicate;
    CHK_THROW_MESSAGE(VecDuplicate(vec, &duplicate), "Failed to duplicate Vec");
    return SmartPetscObj<Vec>(duplicate);
  } else {
    return SmartPetscObj<Vec>();
  }
};

inline SmartPetscObj<Vec> smartVectorDuplicate(Vec &vec) {
  Vec duplicate;
  CHK_THROW_MESSAGE(VecDuplicate(vec, &duplicate), "Failed to duplicate Vec");
  return SmartPetscObj<Vec>(duplicate);
};

inline SmartPetscObj<Mat> smartMatDuplicate(Mat &mat, MatDuplicateOption op) {
  Mat duplicate;
  CHK_THROW_MESSAGE(MatDuplicate(mat, op, &duplicate),
                    "Failed to duplicate Mat");
  return SmartPetscObj<Mat>(duplicate);
};

inline SmartPetscObj<Mat> smartMatDuplicate(SmartPetscObj<Mat> &mat,
                                            MatDuplicateOption op) {
  if (mat.use_count()) {
    Mat duplicate;
    CHK_THROW_MESSAGE(MatDuplicate(mat, op, &duplicate),
                      "Failed to duplicate Mat");
    return SmartPetscObj<Mat>(duplicate);
  } else {
    return SmartPetscObj<Mat>();
  }
};

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

} // namespace MoFEM

#endif