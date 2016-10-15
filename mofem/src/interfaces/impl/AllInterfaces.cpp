const static int debug = 1;

#include "impl/MeshsetsManager.cpp"
#include "impl/CoordSystemsManager.cpp"
#include "impl/MeshRefinement.cpp"
#include "impl/PrismInterface.cpp"
#include "impl/SeriesRecorderCore.cpp"
#include "impl/NodeMerger.cpp"
#include "impl/BitLevelCoupler.cpp"
#include "impl/PrismsFromSurfaceInterface.cpp"

#ifdef WITH_TETGEN
  #include "impl/TetGenInterface.cpp"
#endif

#ifdef WITH_MED
  #include "impl/MedInterface.cpp"
#endif
