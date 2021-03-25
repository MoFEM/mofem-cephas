/** \file AuxPETSc.hpp
 * \brief Auxiliary MoFEM-PETSc structures
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

struct PairNameFEMethodPtr : public std::pair<std::string, FEMethod *> {

  PairNameFEMethodPtr(std::string name, FEMethod *ptr)
      : std::pair<std::string, FEMethod *>(name, ptr) {}
  template <typename FEMETHOD>
  PairNameFEMethodPtr(std::string name, boost::shared_ptr<FEMETHOD> ptr)
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
  template <typename BASICMETHOD>
  BasicMethodPtr(boost::shared_ptr<BASICMETHOD> ptr)
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

typedef std::deque<PairNameFEMethodPtr> FEMethodsSequence;
typedef std::deque<BasicMethodPtr> BasicMethodsSequence;

} // namespace MoFEM

#endif // __AUXPETSC_HPP__

/**
 * \defgroup mofem_petsc_solvers PETSc solvers
 * \brief PETSc solvers
 *
 * \ingroup mofem
 */
