/**
 * \file poisson_3d_cutfem.cpp
 * \example poisson_3d_cutfem.cpp
 *
 * Example of implementation for discontinuous Galerkin.
 * ccache make && ./poisson_cutfem -file_name /home/karol/mofem_install/mofem-cephas/mofem/tutorials/scl-13/meshes/hex_mesh.cub -immersed_mesh elipse2.stl  && vtkmake
 */

#include <MoFEM.hpp>

constexpr int BASE_DIM = 1;
constexpr int FIELD_DIM = 1;
constexpr int SPACE_DIM = 3;

using namespace MoFEM;

using EntData = EntitiesFieldData::EntData;

using DomainEle = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::DomainEle;
using DomainEleOp = DomainEle::UserDataOperator;

using BoundaryEle =
    PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::BoundaryEle;
using BoundaryEleOp = BoundaryEle::UserDataOperator;
using FaceSideEle =
    PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::FaceSideEle;
using FaceSideOp = FaceSideEle::UserDataOperator;

using PostProcEle =  PostProcBrokenMeshInMoab<DomainEle>;

static double penalty = 1e6;
static double phi = -1; // 1 - symmetric Nitsche, 0 - nonsymmetric, -1 antisymmetrica
static double nitsche = 1;
static unsigned int max_ref_level = 5;

#include <EntityRefine.hpp>
#include <PoissonCutFEM.hpp>

using OpDomainGradGrad = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::BiLinearForm<GAUSS>::OpGradGrad<BASE_DIM, FIELD_DIM, SPACE_DIM>;
using OpDomainSource = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpSource<BASE_DIM, FIELD_DIM>;

PetscBool is_test = PETSC_FALSE;

auto u_exact = [](const double x, const double y, const double) {
  if (is_test)
    return x * x * y * y;
  else
    return cos(2 * x * M_PI) * cos(2 * y * M_PI);
};

auto u_grad_exact = [](const double x, const double y, const double) {
  if (is_test)
    return FTensor::Tensor1<double, 2>{2 * x * y * y, 2 * x * x * y};
  else
    return FTensor::Tensor1<double, 2>{

        -2 * M_PI * cos(2 * M_PI * y) * sin(2 * M_PI * x),
        -2 * M_PI * cos(2 * M_PI * x) * sin(2 * M_PI * y)

    };
};

auto source = [](const double x, const double y, const double) {
  if (is_test)
    return -(2 * x * x + 2 * y * y);
  else
    return 8 * M_PI * M_PI * cos(2 * x * M_PI) * cos(2 * y * M_PI);
};

using namespace MoFEM;
using namespace Poisson3DCutFEMOperators;

static char help[] = "...\n\n";

struct Poisson3DCutFEM {
  MoFEM::Core *core;

public:
  Poisson3DCutFEM(MoFEM::Interface &m_field);

  // Declaration of the main function to run analysis
  MoFEMErrorCode runProgram();

private:
  // Declaration of other main functions called in runProgram()
  MoFEMErrorCode readMesh();
  MoFEMErrorCode setupProblem();
  MoFEMErrorCode boundaryCondition();
  MoFEMErrorCode assembleSystem();
  MoFEMErrorCode setIntegrationRules();
  MoFEMErrorCode solveSystem();
  // MoFEMErrorCode checkResults();
  MoFEMErrorCode outputResults();

  // MoFEM interfaces
  MoFEM::Interface &mField;
  Simple *simpleInterface;

  moab::Core coreIm;
  moab::Interface &moabIm;

  // Field name and approximation order
  std::string domainField;
  int oRder;
};

Poisson3DCutFEM::Poisson3DCutFEM(MoFEM::Interface &m_field)
    : domainField("U"), mField(m_field), moabIm(coreIm), oRder(4) {}

//! [Read mesh]
MoFEMErrorCode Poisson3DCutFEM::readMesh() {
  MoFEMFunctionBegin;

  CHKERR mField.getInterface(simpleInterface);
  CHKERR simpleInterface->getOptions();
  CHKERR simpleInterface->loadFile();

  char immersed_mesh_file_name[255];
  CHKERR PetscOptionsGetString(PETSC_NULL, "", "-immersed_mesh", immersed_mesh_file_name,
                               255, PETSC_NULL);
  CHKERR moabIm.load_file(immersed_mesh_file_name, 0, "");

  MoFEMFunctionReturn(0);
}
//! [Read mesh]

//! [Setup problem]
MoFEMErrorCode Poisson3DCutFEM::setupProblem() {
  MoFEMFunctionBegin;

  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &oRder, PETSC_NULL);
  CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-penalty", &penalty,
                               PETSC_NULL);
  CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-phi", &phi, PETSC_NULL);
  CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-nitsche", &nitsche,
                               PETSC_NULL);
  PetscOptionsGetBool(PETSC_NULL, "", "-is_test", &is_test, PETSC_NULL);

  MOFEM_LOG("WORLD", Sev::inform) << "Set order: " << oRder;
  MOFEM_LOG("WORLD", Sev::inform) << "Set penalty: " << penalty;
  MOFEM_LOG("WORLD", Sev::inform) << "Set phi: " << phi;
  MOFEM_LOG("WORLD", Sev::inform) << "Set nitche: " << nitsche;
  MOFEM_LOG("WORLD", Sev::inform) << "Set test: " << (is_test == PETSC_TRUE);

  CHKERR simpleInterface->addDomainField(domainField, H1, DEMKOWICZ_JACOBI_BASE,
                                         1);
  CHKERR simpleInterface->setFieldOrder(domainField, oRder);

  // This is only for debugging and experimentation, to see boundary and edge
  // elements.
  auto save_shared = [&](auto meshset, std::string prefix) {
    MoFEMFunctionBegin;
    auto file_name =
        prefix + "_" +
        boost::lexical_cast<std::string>(mField.get_comm_rank()) + ".vtk";
    CHKERR mField.get_moab().write_file(file_name.c_str(), "VTK", "", &meshset,
                                        1);
    MoFEMFunctionReturn(0);
  };

  // CHKERR save_shared(simpleInterface->getBoundaryMeshSet(), "bdy");
  // CHKERR save_shared(simpleInterface->getSkeletonMeshSet(), "skel");

  CHKERR simpleInterface->addDataField("DISTANCE_FROM_SURFACE", H1,
                                         DEMKOWICZ_JACOBI_BASE, 1);
  CHKERR simpleInterface->setFieldOrder("DISTANCE_FROM_SURFACE", 2);
  CHKERR simpleInterface->setUp();

  Range immersed_mesh_ents;
  SurfaceKDTree surface_distance(mField, moabIm);
  {
    CHKERR moabIm.get_entities_by_dimension(0, 2, immersed_mesh_ents, false);
    std::cout << " We are at a problem at line: " << immersed_mesh_ents << "\n";
    CHKERR surface_distance.takeASkin(immersed_mesh_ents);
    CHKERR surface_distance.buildTree(immersed_mesh_ents);
    std::cout << " We are at a problem at line: 2" << immersed_mesh_ents << "\n";
    CHKERR surface_distance.setDistanceFromSurface("DISTANCE_FROM_SURFACE");
    CHKERR surface_distance.copySurface(immersed_mesh_ents, mField.get_moab());

    auto set_vets = [&](boost::shared_ptr<FieldEntity> ent_ptr) {
      MoFEMFunctionBegin;

      EntityHandle ent = ent_ptr->getEnt();
      double coords[3];
      CHKERR mField.get_moab().get_coords(&ent, 1, coords);
      std::cout << " We are at a problem at line: " << __LINE__ << "\n";
      double delta[3];
      CHKERR surface_distance.findClosestPointToTheSurface(
          coords[0], coords[1], coords[2], delta[0], delta[1], delta[2]);
      ent_ptr->getEntFieldData()[0] = delta[0];
      ent_ptr->getEntFieldData()[1] = delta[1];
      ent_ptr->getEntFieldData()[2] = delta[2];
      MoFEMFunctionReturn(0);
    };
  }

  
  EntityHandle ref_level_meshset = 0;

  CHKERR save_shared(ref_level_meshset, "sdf_check");

  auto cut_mesh = mField.getInterface<CutMeshInterface>();
  CHKERR cut_mesh->setSurface(immersed_mesh_ents);
  // Create tag storing nodal positions
  double def_position[] = {0, 0, 0};
  Tag th;
  CHKERR mField.get_moab().tag_get_handle("POSITION", 3, MB_TYPE_DOUBLE, th,
                             MB_TAG_CREAT | MB_TAG_SPARSE, def_position);
  // Set tag values with coordinates of nodes
  CHKERR cut_mesh->setTagData(th);
  CHKERR cut_mesh->buildTree();

  Range hexes;
  CHKERR mField.get_moab().get_entities_by_dimension(0, 3, hexes, false);

  auto cut_mesh_intersect = [&](Range const &ents, Range &intersect_edges, int cbit) {
    MoFEMFunctionBeginHot;
    CHKERR cut_mesh->setVolume(ents);
    CHKERR cut_mesh->createSurfaceLevelSets();

    CHKERR cut_mesh->findEdgesToCut(ents);
    // CHKERR cut_mesh->cutOnly(ents, cbit, th, 1e-4, 1e-4);
    Range const &cut_edges = cut_mesh->getCutEdges();
    intersect_edges = cut_edges;

    CHKERR CutMeshInterface::SaveData(
        mField.get_moab())("cut_edges_" + std::to_string(cbit) + ".vtk",
        cut_edges);

    Range cut_faces;
    CHKERR mField.get_moab().get_adjacencies(
        cut_edges, 2, true, cut_faces, moab::Interface::UNION);
    CHKERR CutMeshInterface::SaveData(mField.get_moab())(
        "cut_faces_" + std::to_string(cbit) + ".vtk", cut_faces);
    // intersect_edges.merge(cut_edges);
    // std::cout << " intersect_edges " << intersect_vols << "\n";
    MoFEMFunctionReturnHot(0);
  };

  CHKERR CutMeshInterface::SaveData(mField.get_moab())("intersect_vols0.vtk",
                                                       hexes);

  Range intersect_edges, intersect_vols, intersect_faces;
  CHKERR cut_mesh_intersect(hexes, intersect_edges, 1);
  CHKERR mField.get_moab().get_adjacencies(
      intersect_edges, 3, false, intersect_vols, moab::Interface::UNION);
  CHKERR mField.get_moab().get_adjacencies(
      intersect_edges, 2, false, intersect_faces, moab::Interface::UNION);

  auto bit = [](auto l) { return BitRefLevel().set(l); };

  auto refine_tet_level = [&](Range &level_ents, Range const &cut_edges, int const &l) {
    MoFEMFunctionBegin;
    BitRefManager *bit_ref_manager;
    MeshRefinement *refine = mField.getInterface<MeshRefinement>();
    CHKERR mField.getInterface(bit_ref_manager);
    // Range level_ents;
    // CHKERR bit_ref_manager->getEntitiesByDimAndRefLevel(
    //     bit(l - 1), BitRefLevel().set(), 3, level_ents);
    // if (l <= 2) {
    if(false ) {
      // Range edges;
      // CHKERR mField.get_moab().get_adjacencies(tets.subset_by_dimension(3), 1, true, edges,
      //                                          moab::Interface::UNION);

      // EntityHandle meshset_ref_edges;
      // CHKERR mField.get_moab().create_meshset(MESHSET_SET, meshset_ref_edges);
      // CHKERR mField.get_moab().add_entities(meshset_ref_edges, edges);

      // CHKERR refine->addVerticesInTheMiddleOfEdges(meshset_ref_edges, bit(l + 0));

      // CHKERR refine->refineTets(tets, bit(l + 0));
    } else {
      // EntityHandle meshset_ref_edges;
      // CHKERR mField.get_moab().create_meshset(MESHSET_SET, meshset_ref_edges);
      // CHKERR mField.get_moab().add_entities(meshset_ref_edges, cut_edges);
      // CHKERR refine->addVerticesInTheMiddleOfEdges(meshset_ref_edges, bit(l + 0));

      // if (l <= 2)
      //   CHKERR refine->refineTets(level_ents, bit(l + 0));
      // else 
      {
        Range tets_to_refine, all_edges;
        CHKERR mField.get_moab().get_adjacencies(cut_edges, 3, false, tets_to_refine,
                                                 moab::Interface::UNION);
        CHKERR mField.get_moab().get_adjacencies(tets_to_refine, 1, true, all_edges,
                                                 moab::Interface::UNION);
        EntityHandle meshset_ref_edges;
        CHKERR mField.get_moab().create_meshset(MESHSET_SET, meshset_ref_edges);
        CHKERR mField.get_moab().add_entities(meshset_ref_edges, all_edges);
  
        // CHKERR bit_ref_manager->setNthBitRefLevel(level_ents, l - 1, l);
        // Range intersect_vols = intersect(tets_to_refine, level_ents);
        // CHKERR bit_ref_manager->setEntitiesBitRefLevel(level_ents, bit(l - 1));
        
        //FIXME:
        level_ents = subtract(level_ents, tets_to_refine);
        // CHKERR bit_ref_manager->setEntitiesBitRefLevel(level_ents, bit(l - 1));
        // level_ents = tets_to_refine;
        CHKERR bit_ref_manager->setEntitiesBitRefLevel(level_ents, bit(l));

        CHKERR refine->addVerticesInTheMiddleOfEdges(
            meshset_ref_edges, bit(l + 0));

        CHKERR refine->refineTetsHangingNodes(tets_to_refine, bit(l));
      }
    }

    MoFEMFunctionReturn(0);
  };

  auto bit_mng = mField.getInterface<BitRefManager>();
  CHKERR bit_mng->setBitRefLevelByDim(0, 3, bit(0));
  HexMeshRefinement hex_refine(*core);
  CHKERR hex_refine.refineHexToTets(intersect_vols, bit(0 + 1));
  constexpr bool debug = false;

  // bit 0 - hexes
  // bit 1 - tets
  // bit 2 - refined tets
  // bit 3 - more refined tets ...

  Range level0;
  CHKERR bit_mng->getEntitiesByDimAndRefLevel(bit(0), BitRefLevel().set(), 3,
                                              level0);
  level0 = subtract(level0, intersect_vols);
  CHKERR bit_mng->setEntitiesBitRefLevel(level0, bit(1));
  // CHKERR bit_mng->setEntitiesBitRefLevel(intersect_vols, bit(1));
  // CHKERR CutMeshInterface::SaveData(mField.get_moab())("intersect_vols0.vtk",
  //                                                      level0);

  Range level1;
  CHKERR bit_mng->getEntitiesByDimAndRefLevel(bit(1), BitRefLevel().set(), 3, level1);
  CHKERR CutMeshInterface::SaveData(mField.get_moab())("intersect_vols1.vtk",
                                                       level1);
  auto refine_cut_tets = [&](int bits) {
    MoFEMFunctionBegin;

    for (int l = 2; l < bits + 1; l++)
    {
      Range level_ents;
      CHKERR bit_mng->getEntitiesByDimAndRefLevel(bit(l - 1), BitRefLevel().set(), 3, level_ents);
      Range cut_edges;
      CHKERR cut_mesh_intersect(level_ents, cut_edges, l);
      CHKERR refine_tet_level(level_ents, cut_edges, l);
      Range refined_tets;
      CHKERR bit_mng->getEntitiesByDimAndRefLevel(bit(l), BitRefLevel().set(), 3, refined_tets);
      CHKERR CutMeshInterface::SaveData(mField.get_moab())("intersect_vols" + std::to_string(l) + ".vtk", refined_tets);
    }
    
    MoFEMFunctionReturn(0);
  };

  CHKERR refine_cut_tets(max_ref_level);

  std::function<Range(Range &)> get_all_children = [&](Range &ents) -> Range {
    Range all_children;
    CHKERR(bit_mng->updateRangeByChildren(ents, all_children));
    // if (!all_children.empty()) {
    //   ents.merge(get_all_children(all_children));
    // }

    Range children_level;
    for (auto &ent : all_children) {
      Range parent(ent, ent);
      children_level.merge(get_all_children(parent));
    }

    return children_level.empty() ? ents : children_level;
  };

  std::function<Range(const EntityHandle &)> get_all_children_from_ent =
      [&](const EntityHandle &ent) -> Range {
    Range all_children(ent, ent);
    CHKERR(bit_mng->updateRangeByChildren(all_children, all_children));
    for (auto &enit : all_children)
      all_children.merge(get_all_children_from_ent(enit));

    return all_children.empty() ? Range(ent, ent) : all_children;
  };

  for (auto &ent : intersect_vols) {
    Range single_ent(ent, ent);
    Range all_refined_ents = get_all_children(single_ent);
    all_refined_ents = all_refined_ents.subset_by_dimension(3);
    // for(int = 0; i < max_ref_level; i++)
    // {
    //   CHKERR bit_mng->updateRangeByChildren(single_ent, single_ent);
    // }
    CHKERR CutMeshInterface::SaveData(mField.get_moab())(
        "refined_tet" + std::to_string(1) + ".vtk", all_refined_ents);
    
    CHKERR bit_mng->filterEntitiesByRefLevel(
        bit(max_ref_level), BitRefLevel().set(), all_refined_ents);

    CHKERR CutMeshInterface::SaveData(mField.get_moab())(
        "refined_tet" + std::to_string(2) + ".vtk", all_refined_ents);
    break;
  }

  // CHKERR CutMeshInterface::SaveData(mField.get_moab())(
  //     "refined_tet" + std::to_string(3) + ".vtk", all_refined_ents);
  MoFEMFunctionReturn(0);
}
//! [Setup problem]

//! [Boundary condition]
MoFEMErrorCode Poisson3DCutFEM::boundaryCondition() {
  MoFEMFunctionBegin;
  //FIXME: Implement boundary condition
  MoFEMFunctionReturn(0);
}

//! [Boundary condition]

//! [Assemble system]
MoFEMErrorCode Poisson3DCutFEM::assembleSystem() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();

  auto add_base_ops = [&](auto &pipeline) {
    auto det_ptr = boost::make_shared<VectorDouble>();
    auto jac_ptr = boost::make_shared<MatrixDouble>();
    auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
    pipeline.push_back(new OpCalculateHOJacForFace(jac_ptr));
    pipeline.push_back(new OpInvertMatrix<2>(jac_ptr, det_ptr, inv_jac_ptr));
    pipeline.push_back(new OpSetInvJacL2ForFace(inv_jac_ptr));
  };

  add_base_ops(pipeline_mng->getOpDomainLhsPipeline());
  pipeline_mng->getOpDomainLhsPipeline().push_back(new OpDomainGradGrad(
      domainField, domainField,
      [](const double, const double, const double) { return 1; }));
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpDomainSource(domainField, source));

  MoFEMFunctionReturn(0);
}
//! [Assemble system]

//! [Set integration rules]
MoFEMErrorCode Poisson3DCutFEM::setIntegrationRules() {
  MoFEMFunctionBegin;

  auto rule_lhs = [](int, int, int p) -> int { return 2 * p; };
  auto rule_rhs = [](int, int, int p) -> int { return 2 * p; };
  // auto rule_2 = [this](int, int, int) { return 2 * oRder; };

  auto pipeline_mng = mField.getInterface<PipelineManager>();
  CHKERR pipeline_mng->setDomainLhsIntegrationRule(rule_lhs);
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule_rhs);

  // CHKERR pipeline_mng->setSkeletonLhsIntegrationRule(rule_2);
  // CHKERR pipeline_mng->setSkeletonRhsIntegrationRule(rule_2);
  // CHKERR pipeline_mng->setBoundaryLhsIntegrationRule(rule_2);
  // CHKERR pipeline_mng->setBoundaryRhsIntegrationRule(rule_2);

  // pipeline_mng->getDomainRhsFE()->getRuleHook = [](int, int, int) { return -1; };
  // pipeline_mng->getDomainRhsFE()->setRuleHook = [](int, int, int) { return -1; };
  // pipeline_mng->getSkeletonRhsFE()->getRuleHook = [](int, int, int) { return -1; };
  // pipeline_mng->getSkeletonRhsFE()->setRuleHook = [](int, int, int) { return -1; };
  // fe->getRuleHook = [](int, int, int) { return -1; }; // VolRule();
  // fe->setRuleHook = SetIntegrationAtFrontVolume(frontVertices,
  // frontAdjEdges);

  MoFEMFunctionReturn(0);
}
//! [Set integration rules]

//! [Solve system]
MoFEMErrorCode Poisson3DCutFEM::solveSystem() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();

  auto ksp_solver = pipeline_mng->createKSP();
  CHKERR KSPSetFromOptions(ksp_solver);
  // CHKERR KSPSetUp(ksp_solver);

  // Create RHS and solution vectors
  auto dm = simpleInterface->getDM();
  auto F = createDMVector(dm);
  auto D = vectorDuplicate(F);

  // CHKERR KSPSolve(ksp_solver, F, D);

  // Scatter result data on the mesh
  CHKERR VecGhostUpdateBegin(D, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(D, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR DMoFEMMeshToLocalVector(dm, D, INSERT_VALUES, SCATTER_REVERSE);

  MoFEMFunctionReturn(0);
}
//! [Solve system]

//! [Output results]
MoFEMErrorCode Poisson3DCutFEM::outputResults() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();
  pipeline_mng->getDomainLhsFE().reset();
  pipeline_mng->getSkeletonRhsFE().reset();
  pipeline_mng->getSkeletonLhsFE().reset();
  pipeline_mng->getBoundaryRhsFE().reset();
  pipeline_mng->getBoundaryLhsFE().reset();

  auto post_proc_fe = boost::make_shared<PostProcEle>(mField);

  auto u_ptr = boost::make_shared<VectorDouble>();
  post_proc_fe->getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues(domainField, u_ptr));

  // auto md_ptr = boost::make_shared<MatrixDouble>();
  auto md_ptr = boost::make_shared<VectorDouble>();
  post_proc_fe->getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues("DISTANCE_FROM_SURFACE", md_ptr));

  using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;
  std::cout << "We are at a POSTPROCESSING at line: " << __LINE__ << "\n";

  auto grad_ptr = boost::make_shared<MatrixDouble>();
  post_proc_fe->getOpPtrVector().push_back(new OpCalculateScalarFieldGradient<3>("DISTANCE_FROM_SURFACE", grad_ptr));

  post_proc_fe->getOpPtrVector().push_back(

      new OpPPMap(

          post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

          {{"U", u_ptr}, {"DISTANCE_FROM_SURFACE", md_ptr}},

          {{"DISTANCE_FROM_SURFACE_GRAD", grad_ptr}},

          {},

          {})

  );

  pipeline_mng->getDomainRhsFE() = post_proc_fe;
  CHKERR pipeline_mng->loopFiniteElements();
  CHKERR post_proc_fe->writeFile("out_result.h5m");

  MoFEMFunctionReturn(0);
}
//! [Output results]

//! [Run program]
MoFEMErrorCode Poisson3DCutFEM::runProgram() {
  MoFEMFunctionBegin;

  CHKERR readMesh();
  CHKERR setupProblem();
  CHKERR boundaryCondition();
  CHKERR assembleSystem();
  CHKERR setIntegrationRules();
  CHKERR solveSystem();
  CHKERR outputResults();
  // CHKERR checkResults();

  MoFEMFunctionReturn(0);
}
//! [Run program]

//! [Main]
int main(int argc, char *argv[]) {

  // Initialisation of MoFEM/PETSc and MOAB data structures
  const char param_file[] = "param_file.petsc";
  MoFEM::Core::Initialize(&argc, &argv, param_file, help);

  // Error handling
  try {
    // Register MoFEM discrete manager in PETSc
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Create MOAB instance
    moab::Core mb_instance;              // mesh database
    moab::Interface &moab = mb_instance; // mesh database interface

    // Create MoFEM instance
    MoFEM::Core core(moab);           // finite element database
    MoFEM::Interface &m_field = core; // finite element interface

    // Run the main analysis
    Poisson3DCutFEM poisson_problem(m_field);
    poisson_problem.core = &core;
    CHKERR poisson_problem.runProgram();
  }
  CATCH_ERRORS;

  // Finish work: cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}
//! [Main]