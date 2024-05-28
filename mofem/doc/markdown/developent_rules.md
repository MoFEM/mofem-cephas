# Development rules

- MoFEM core has to allow for evolution on the future. Need to be modular, and extendable, remain small as possible. Free from implementation of specific physics.
- MoFEM Users modules should present good development practices. Examples of
implementation for particular physical problems.
- User module is specific code, with bespoke solutions, to be used by external 
users. Solutions can be close, fit to problem.


## General rules

- If code part do not have be fast, make it generic.
- Make generic code as simple as possible.
- Do not duplicate code with similar functionality.
- Make code fast. To check what need to be fast profile code.

## Hide complexities.

- Isolate complexities. You can focus only on one thing. 
- Avoid bool flags and code branching. That makes code default to debug and understand.
- Perfect operator is a black box. 
- Black box do not create dependencies. Enable future evolution and development.
- Do not expose members of user data operators. Make members of the function private. 
- Common data for series of operators are local. 
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

## Avoid decencies

- Do not add software dependencies to core library.
- If you have to add dependency, make it optional.
- More libraries bigger likelihood that code will break at compilation.

