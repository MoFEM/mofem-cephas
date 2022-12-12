const static int debug = 0;

#include <MoFEM.hpp>

#include "impl/ProblemsManager.cpp"
#include "impl/Simple.cpp"
#include "impl/OperatorsTester.cpp"
#include "impl/PipelineManager.cpp"
#include "impl/MeshsetsManager.cpp"
#include "impl/MeshRefinement.cpp"
#include "impl/PrismInterface.cpp"
#include "impl/SeriesRecorder.cpp"
#include "impl/NodeMerger.cpp"
#include "impl/PrismsFromSurfaceInterface.cpp"
#include "impl/FieldEvaluator.cpp"
#include "impl/Tools.cpp"
#include "impl/BcManager.cpp"

#ifdef WITH_TETGEN
  #include "impl/TetGenInterface.cpp"
#endif

#ifdef WITH_MED
  #include "impl/MedInterface.cpp"
#endif

#include "impl/CutMeshInterface.cpp"
