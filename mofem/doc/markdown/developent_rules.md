# Development rules

- MoFEM core is to be modular, and extendable, remain small as possible. That will enable its evolution in future.
- MoFEM core is free from implementation of specific physics. 	
- MoFEM Users modules should present good development practices. 
- MoFEM Users modules are examples of implementation for particular physical problems.
- User module is specific code, with bespoke solutions, to be used by external 
users. Solutions can be close, bespoke, tailored to problem.
- Code always works in parallel.

## General rules

- If code part do not have be fast, make it generic.
- Make generic code as simple as possible.
- Do not duplicate code with similar functionality.
- Make code fast. To check what need to be fast profile code.

## Hide complexities.

- Isolate complexities. You can focus only on one thing. 
- Avoid bool flags and code branching. That makes code difficult to test, debug and understand.
- Perfect operator is a black box. 
- Black box do not create dependencies. Enable future evolution and development.
- Do not expose members of user data operators. Make members of the function private. 
- Common data for series of operators are local. 

```c++
template <int DIM, AssemblyType A, IntegrationType I, typename BoundaryEleOp>
MoFEMErrorCode opFactoryBoundaryRhs(
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string sigma, std::string u, bool is_axisymmetric = false) {
  MoFEMFunctionBegin;

  using C = ContactIntegrators<BoundaryEleOp>;

  auto common_data_ptr = boost::make_shared<ContactOps::CommonData>();

  pip.push_back(new OpCalculateVectorFieldValues<DIM>(
      u, common_data_ptr->contactDispPtr()));
  pip.push_back(new OpCalculateHVecTensorTrace<DIM, BoundaryEleOp>(
      sigma, common_data_ptr->contactTractionPtr()));
  pip.push_back(
      new typename C::template Assembly<A>::template OpConstrainBoundaryRhs<
          DIM, GAUSS>(sigma, common_data_ptr, is_axisymmetric));

  MoFEMFunctionReturn(0);
}
```

- Prefect operator is dimension, and element type agnostic.
- Push operators in factory function, or lambda function if ad hoc solution.
- Push operator by operator defined in "cpp" file. Hide implementation.

hpp file:
```c++
namespace BlackBox  {
struct OpBlackBox; // only declaration
OpBlackBox* createOpBlackBox();
}
```

cpp: file:
```c++
namespace BlackBox  {
struct OpBalckBox: private ForcesAndSourcesCore::UserDataOperator {
	private:
	// definition
};
OpBlackBox* createOpBlackBox() { return new OpBalckBox(); }
}
```

## Testing

- Test code not method.
- Tests are fast.

## Avoid dependencies

- Do not add software dependencies to core library.
- If you have to add dependency, make it optional.
- More libraries bigger likelihood that code will break at compilation.
- Some decencies should be restricted to users modules
- All dependencies are managed in Spack

