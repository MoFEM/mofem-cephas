# Development rules

- MoFEM core is to be modular, extendable, and to remain as small as possible to enable its future evolution.
- MoFEM core is free from implementation of specific physics. 	
- MoFEM User modules should present good development practices. 
- MoFEM User modules are examples of implementation for particular physical problems.
- User modules is specific code, with bespoke solutions, to be used by external 
users. Solutions can be closed, bespoke, and/or tailored.
- Code always works in parallel.

## General rules

- If a section of the code does not have to be fast, make it generic.
- Make generic code as simple as possible.
- Do not duplicate code with similar functionality.
- Make code fast. Profile the code to check what needs to be made fast.


## Hide complexities.

- Isolate complexities. You can focus only on one thing. 
- Avoid bool flags and code branching. This makes code difficult to test, debug, and understand.
- A perfect operator is a black box.
- Black boxes do not create dependencies, enabling future evolution and development
- Do not expose members of user data operators. Make members of the function private. 
- Common data for a series of operators are local. 

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

- The perfect operator is dimension and element type agnostic.
- Push operators in a factory function or a lambda function if an ad hoc solution
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
struct OpBlackBox: private ForcesAndSourcesCore::UserDataOperator {
	private:
	// definition
};
OpBlackBox* createOpBlackBox() { return new OpBalckBox(); }
}
```

## Testing

- Test the code, not the method.
- Tests are fast.

## Avoid dependencies

- Do not add software dependencies to core library.
- If you have to add a dependency, make it optional.
- With more libraries, there is a larger likelihood of the code breaking at compilation.
- Some dependencies should be restricted to users modules
- All dependencies are managed managed Spack

