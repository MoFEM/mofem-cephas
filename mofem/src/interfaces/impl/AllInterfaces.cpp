const static int debug = 1;

#include <version.h>
#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <FTensor.hpp>
#include <fem_tools.h>
#include <h1_hdiv_hcurl_l2.h>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <UnknownInterface.hpp>
#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <Core.hpp>

// Interfaces
#include <ProblemsManager.hpp>
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
#include <UpdateMeshsetsAndRanges.hpp>
#include <Tools.hpp>

#include <boost/scoped_ptr.hpp>
#include <moab/AdaptiveKDTree.hpp>
#include <BitLevelCoupler.hpp>

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
