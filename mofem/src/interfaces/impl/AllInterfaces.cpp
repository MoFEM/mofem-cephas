const static int debug = 0;

#include <MoFEM.hpp>

// Interfaces
#include <ProblemsManager.hpp>
#include <MatrixManager.hpp>
#include <Simple.hpp>
#include <ISManager.hpp>
#include <BitRefManager.hpp>
#include <VecManager.hpp>
#include <FieldBlas.hpp>
#include <MeshRefinement.hpp>
#include <SeriesRecorder.hpp>
#include <PrismInterface.hpp>
#include <CutMeshInterface.hpp>
#include <MeshsetsManager.hpp>
#include <CoordSystemsManager.hpp>
#include <TetGenInterface.hpp>
#include <MedInterface.hpp>
#include <NodeMerger.hpp>
#include <PrismsFromSurfaceInterface.hpp>
#include <Tools.hpp>
#include <CommInterface.hpp>

#include "impl/ProblemsManager.cpp"
#include "impl/Simple.cpp"
#include "impl/MeshsetsManager.cpp"
#include "impl/CoordSystemsManager.cpp"
#include "impl/MeshRefinement.cpp"
#include "impl/PrismInterface.cpp"
#include "impl/SeriesRecorder.cpp"
#include "impl/NodeMerger.cpp"
#include "impl/BitLevelCoupler.cpp"
#include "impl/PrismsFromSurfaceInterface.cpp"

#ifdef WITH_TETGEN
  #include "impl/TetGenInterface.cpp"
#endif

#ifdef WITH_MED
  #include "impl/MedInterface.cpp"
#endif

#include "impl/CutMeshInterface.cpp"
