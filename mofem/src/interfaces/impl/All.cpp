// Trick to make it compile faster

const static int debug = 1;

#include <MoFEM.hpp>

#include "impl/ErrorHandler.cpp"
#include "impl/LoopMethods.cpp"
#include "impl/Core.cpp"
#include "impl/DeprecatedCoreInterface.cpp"
#include "impl/FieldCore.cpp"
#include "impl/FECore.cpp"
#include "impl/ProblemsCore.cpp"
#include "impl/MatrixManager.cpp"
#include "impl/CommCore.cpp"
#include "impl/DeleteCore.cpp"
#include "impl/ISManager.cpp"
#include "impl/VecManager.cpp"
#include "impl/FieldBlas.cpp"
#include "impl/BitRefManager.cpp"
#include "impl/Tools.cpp"
#include "impl/FieldEvaluator.cpp"
