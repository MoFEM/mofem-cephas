/** \file AuxPETSc.hpp
 * \brief Auxiliary MoFEM-PETSc structures
 */



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
