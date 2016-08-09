const static int debug = 1;

#include "impl/MeshRefinment.cpp"
#include "impl/PrismInterfaceCore.cpp"
#include "impl/SeriesRecorderCore.cpp"
#include "impl/NodeMerger.cpp"
#include "impl/BitLevelCoupler.cpp"
#include "impl/PrismsFromSurfaceInterface.cpp"

#ifdef WITH_TETGEN
  #include "impl/TetGenInterface.cpp"
#endif

#ifdef WITH_NETGEN
  #include "impl/NetGenInterface.cpp"
#endif
