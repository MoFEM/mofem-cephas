// Trick to make it compile faster

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
#include <Tools.hpp>

#include <boost/scoped_ptr.hpp>
#include <moab/AdaptiveKDTree.hpp>
#include <BitLevelCoupler.hpp>

#include "impl/ErrorHandler.cpp"
#include "impl/LoopMethods.cpp"
#include "impl/Core.cpp"
#include "impl/DeprecatedCoreInterface.cpp"
#include "impl/CreateMat.cpp"
#include "impl/FieldCore.cpp"
#include "impl/FECore.cpp"
#include "impl/ProblemsCore.cpp"
#include "impl/CommCore.cpp"
#include "impl/DeleteCore.cpp"
#include "impl/ISManager.cpp"
#include "impl/VecManager.cpp"
#include "impl/FieldBlas.cpp"
#include "impl/BitRefManager.cpp"
#include "impl/Tools.cpp"
