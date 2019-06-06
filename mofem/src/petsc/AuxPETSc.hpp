/** \file AuxPETSc.hpp
 * \brief Auxuliary MoFEM-PETSc structures
 */

/* This file is part of MoFEM.
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#ifndef __AUXPETSC_HPP__
#define __AUXPETSC_HPP__

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
  PetscErrorCode ierr;
  int cnt;
  ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(obj), &cnt);
  CHKERRABORT(PetscObjectComm(MoFEM::getPetscObject(obj)), ierr);
  if (cnt > 1) {
    ierr = PetscObjectDereference(MoFEM::getPetscObject(obj));
    CHKERRABORT(PetscObjectComm(MoFEM::getPetscObject(obj)), ierr);
  } else {
    ierr = PetscObjectDestroy(reinterpret_cast<PetscObject *>(&obj));
    CHKERRABORT(PetscObjectComm(MoFEM::getPetscObject(obj)), ierr);
  }
}

namespace MoFEM {

struct PairNameFEMethodPtr : public std::pair<std::string, FEMethod *> {

  PairNameFEMethodPtr(std::string name, FEMethod *ptr)
      : std::pair<std::string, FEMethod *>(name, ptr) {}
  PairNameFEMethodPtr(std::string name, boost::shared_ptr<FEMethod> &ptr)
      : std::pair<std::string, FEMethod *>(name, ptr.get()), fePtr(ptr) {}
  virtual ~PairNameFEMethodPtr() {}

  inline boost::shared_ptr<BasicMethod> getSharedPtr() const {
    if (!fePtr)
      THROW_MESSAGE("Shared pointer not set. You have to be using raw "
                    "pointer, that is unsafe.");
    return fePtr;
  }

private:
  boost::shared_ptr<FEMethod> fePtr;
};

struct BasicMethodPtr {
  BasicMethodPtr(BasicMethod *ptr) : rawPtr(ptr) {}
  BasicMethodPtr(boost::shared_ptr<BasicMethod> &ptr)
      : rawPtr(ptr.get()), bmPtr(ptr) {}
  BasicMethodPtr(boost::shared_ptr<FEMethod> &ptr)
      : rawPtr(ptr.get()), bmPtr(ptr) {}
  inline BasicMethod &operator*() const { return *rawPtr; };
  inline BasicMethod *operator->() const { return rawPtr; }

  inline boost::shared_ptr<BasicMethod> getSharedPtr() const {
    if (!bmPtr)
      THROW_MESSAGE("Shared pointer not set. You have to be using raw "
                    "pointer, that is unsafe.");
    return bmPtr;
  }

private:
  BasicMethod *rawPtr;
  boost::shared_ptr<BasicMethod> bmPtr;
};

typedef std::vector<PairNameFEMethodPtr> FEMethodsSequence;
typedef std::vector<BasicMethodPtr> BasicMethodsSequence;

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

  SmartPetscObj()
      : boost::intrusive_ptr<typename std::remove_pointer<OBJ>::type>() {}
  SmartPetscObj(OBJ o)
      : boost::intrusive_ptr<typename std::remove_pointer<OBJ>::type>(o,
                                                                      false) {}
  operator OBJ() { return this->get(); }

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
 * @brief Create duplicate vector of smart vector
 *
 * For details abut arguments see here:
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecDuplicate.html>VecDuplicate</a>.
 */
auto smartVectorDuplicate = [](SmartPetscObj<Vec> &vec) {
  if (vec.use_count()) {
    Vec duplicate;
    ierr = VecDuplicate(vec, &duplicate);
    CHKERRABORT(PETSC_COMM_SELF, ierr);
    return SmartPetscObj<Vec>(duplicate);
  } else {
    return SmartPetscObj<Vec>();
  }
};

} // namespace MoFEM

#endif // __AUXPETSC_HPP__
