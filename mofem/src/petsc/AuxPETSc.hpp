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

// intrusive_ptr has to be global
template <typename OBJ> void intrusive_ptr_add_ref(OBJ obj) {
  CHKERR PetscObjectReference(reinterpret_cast<PetscObject>(obj));
}

template <typename OBJ> void intrusive_ptr_release(OBJ obj) {
  int cnt;
  PetscObjectGetReference(reinterpret_cast<PetscObject>(obj), &cnt);
  if (cnt > 1)
    CHKERR PetscObjectDereference(reinterpret_cast<PetscObject>(obj));
  else
    CHKERR PetscObjectDestroy(reinterpret_cast<PetscObject *>(&obj));
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
      PetscObjectGetReference(reinterpret_cast<PetscObject>(this->get()), &cnt);
      return cnt;
    } else
      return 0;
  }
};

auto get_mofem_dm = [](MPI_Comm comm,
                       const std::string dm_type_name = "MOFEM") {
  DM dm;
  CHKERR DMCreate(comm, &dm);
  CHKERR DMSetType(dm, dm_type_name.c_str());
  return SmartPetscObj<DM>(dm);
};

} // namespace MoFEM

#endif // __AUXPETSC_HPP__
