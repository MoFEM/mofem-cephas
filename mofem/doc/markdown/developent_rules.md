# Development rules

## General rules

- If code part do not have be fast, make it generic.
- Make generic code as simple as possible.
- Do not delicate code with similar funionality
- Make code fast. To check what need to be fast profile code.

## Hide complexities.

- Perfect operator is a black box. Do not expose members of user data operators. Make members of the function private. 
- Common data for series of operators are local. Push operators in factory function, or lambda function if ad hoc solution.
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
struct OpBalckBox: public  ForcesAndSourcesCore::UserDataOperator {
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

