// Trick to make it compile faster

const static int debug = 1;

#include "impl/LoopMethods.cpp"
#include "impl/Core.cpp"
#include "impl/CreateMat.cpp"
#include "impl/Vectors.cpp"
#include "impl/FieldBlas.cpp"
#include "impl/GetAdjacancies.cpp"
#include "impl/FieldCore.cpp"
#include "impl/FECore.cpp"
#include "impl/ProblemsCore.cpp"
#include "impl/CommCore.cpp"
#include "impl/DeleteCore.cpp"
#include "impl/MeshRefinment.cpp"
#include "impl/PrismInterfaceCore.cpp"
#include "impl/SeriesRecorderCore.cpp"
#include "impl/NodeMerger.cpp"
#include "impl/BitLevelCoupler.cpp"
#include "impl/PrismsFromSurfaceInterface.cpp"
