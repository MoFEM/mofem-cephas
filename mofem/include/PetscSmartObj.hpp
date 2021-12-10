/** \file PetscSmartObj.hpp
 * \brief Petsc smart obj declarations
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
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
auto createSmartDM = [](MPI_Comm comm, const std::string dm_type_name) {
  DM dm;
  ierr = DMCreate(comm, &dm);
  CHKERRABORT(comm, ierr);
  ierr = DMSetType(dm, dm_type_name.c_str());
  CHKERRABORT(comm, ierr);
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
  ierr = PetscObjectGetComm(obj, &comm);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
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
auto createSmartGhostVector = [](MPI_Comm comm, PetscInt n, PetscInt N,
                                 PetscInt nghost, const PetscInt ghosts[]) {
  Vec vv;
  ierr = VecCreateGhost(comm, n, N, nghost, ghosts, &vv);
  CHKERRABORT(comm, ierr);
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
auto createSmartVectorMPI = [](MPI_Comm comm, PetscInt n, PetscInt N) {
  Vec vv;
  ierr = VecCreateMPI(comm, n, N, &vv);
  CHKERRABORT(comm, ierr);
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
    ierr = VecDuplicate(vec, &duplicate);
    CHKERRABORT(PETSC_COMM_SELF, ierr);
    return SmartPetscObj<Vec>(duplicate);
  } else {
    return SmartPetscObj<Vec>();
  }
};

inline SmartPetscObj<Vec> smartVectorDuplicate(Vec &vec) {
  Vec duplicate;
  ierr = VecDuplicate(vec, &duplicate);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  return SmartPetscObj<Vec>(duplicate);
};

inline SmartPetscObj<Mat> smartMatDuplicate(Mat &mat, MatDuplicateOption op) {
  Mat duplicate;
  ierr = MatDuplicate(mat, op, &duplicate);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  return SmartPetscObj<Mat>(duplicate);
};

inline SmartPetscObj<Mat> smartMatDuplicate(SmartPetscObj<Mat> &mat,
                                            MatDuplicateOption op) {
  if (mat.use_count()) {
    Mat duplicate;
    ierr = MatDuplicate(mat, op, &duplicate);
    CHKERRABORT(PETSC_COMM_SELF, ierr);
    return SmartPetscObj<Mat>(duplicate);
  } else {
    return SmartPetscObj<Mat>();
  }
};

auto createTS = [](MPI_Comm comm) {
  TS ts;
  ierr = TSCreate(comm, &ts);
  CHKERRABORT(comm, ierr);
  return SmartPetscObj<TS>(ts);
};

auto createSNES = [](MPI_Comm comm) {
  SNES snes;
  ierr = SNESCreate(comm, &snes);
  CHKERRABORT(comm, ierr);
  return SmartPetscObj<SNES>(snes);
};

auto createKSP = [](MPI_Comm comm) {
  KSP ksp;
  ierr = KSPCreate(comm, &ksp);
  CHKERRABORT(comm, ierr);
  return SmartPetscObj<KSP>(ksp);
};

auto createPC = [](MPI_Comm comm) {
  PC pc;
  ierr = PCCreate(comm, &pc);
  CHKERRABORT(comm, ierr);
  return SmartPetscObj<PC>(pc);
};

} // namespace MoFEM

#endif