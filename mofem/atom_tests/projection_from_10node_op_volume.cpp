

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

struct OpVolumeCalculation
    : public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

  Vec volumeVec;
  OpVolumeCalculation(Vec volume_vec)
      : MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(
            "MESH_NODE_POSITIONS", UserDataOperator::OPROW),
        volumeVec(volume_vec) {}

  MoFEMErrorCode doWork(int row_side, EntityType row_type,
                        EntitiesFieldData::EntData &row_data) {
    MoFEMFunctionBegin;

    if (row_type != MBVERTEX)
      MoFEMFunctionReturnHot(0);
    double vol = 0;
    auto t_w = getFTensor0IntegrationWeight();

    int nb_gauss_pts = row_data.getN().size1();
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      vol += getVolume() * t_w;
      ++t_w;
    }

    CHKERR VecSetValue(volumeVec, 0, vol, ADD_VALUES);

    MoFEMFunctionReturn(0);
  }
};

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
#else
    CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-my_file",
                                 mesh_file_name, 255, &flg);
#endif
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    }

    const char *option;
    option = ""; //"PARALLEL=BCAST;";//;DEBUG_IO";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    // auto moab_comm_wrap =
    //     boost::make_shared<WrapMPIComm>(PETSC_COMM_WORLD, false);
    // if (pcomm == NULL)
    //   pcomm = new ParallelComm(&moab, moab_comm_wrap->get_comm());

    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // add filds
    CHKERR m_field.add_field("MESH_NODE_POSITIONS", H1, AINSWORTH_LEGENDRE_BASE,
                             3);

    // add finite elements
    CHKERR m_field.add_finite_element("TET_ELEM");
    CHKERR m_field.modify_finite_element_add_field_row("TET_ELEM",
                                                       "MESH_NODE_POSITIONS");
    CHKERR m_field.modify_finite_element_add_field_col("TET_ELEM",
                                                       "MESH_NODE_POSITIONS");
    CHKERR m_field.modify_finite_element_add_field_data("TET_ELEM",
                                                        "MESH_NODE_POSITIONS");

    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);
    Range tets;
    CHKERR moab.get_entities_by_type(0, MBTET, tets, true);
    Range edges;
    CHKERR moab.get_entities_by_type(0, MBEDGE, edges, true);
    CHKERR m_field.getInterface<BitRefManager>()->setElementsBitRefLevel(edges);

    // add ents to field and set app. order
    CHKERR m_field.add_ents_to_field_by_type(0, MBTET, "MESH_NODE_POSITIONS");
    CHKERR m_field.set_field_order(0, MBVERTEX, "MESH_NODE_POSITIONS", 1);
    CHKERR m_field.set_field_order(0, MBEDGE, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(0, MBTRI, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(0, MBTET, "MESH_NODE_POSITIONS", 2);

    // add finite elements entities
    CHKERR m_field.add_ents_to_finite_element_by_type(tets, MBTET, "TET_ELEM");

    // build fields
    CHKERR m_field.build_fields();
    // build finite elements
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);

    Projection10NodeCoordsOnField ent_method(m_field, "MESH_NODE_POSITIONS");
    CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method);

    CHKERR DMRegister_MoFEM("DMMOFEM");
    auto dM = createDM(m_field.get_comm(), "DMMOFEM");
    CHKERR DMMoFEMCreateMoFEM(dM, &m_field, "TET_PROBLEM", bit_level0);
    CHKERR DMMoFEMAddElement(dM, "TET_ELEM");
    CHKERR DMMoFEMSetIsPartitioned(dM, PETSC_FALSE);
    CHKERR DMSetUp_MoFEM(dM);

    auto vol_vec = createVectorMPI(m_field.get_comm(),
                                        m_field.get_comm_rank() ? 0 : 1, 1);

    auto fe_ptr =
        boost::make_shared<VolumeElementForcesAndSourcesCore>(m_field);
    auto material_grad_mat = boost::make_shared<MatrixDouble>();
    auto material_det_vec = boost::make_shared<VectorDouble>();
    fe_ptr->meshPositionsFieldName = "none";
    fe_ptr->getOpPtrVector().push_back(new OpCalculateVectorFieldGradient<3, 3>(
        "MESH_NODE_POSITIONS", material_grad_mat));
    fe_ptr->getOpPtrVector().push_back(
        new OpInvertMatrix<3>(material_grad_mat, material_det_vec, nullptr));
    fe_ptr->getOpPtrVector().push_back(new OpSetHOWeights(material_det_vec));
    fe_ptr->getOpPtrVector().push_back(new OpVolumeCalculation(vol_vec));

    CHKERR DMoFEMLoopFiniteElements(dM, "TET_ELEM", fe_ptr);

    double vol;
    CHKERR VecAssemblyBegin(vol_vec);
    CHKERR VecAssemblyEnd(vol_vec);
    CHKERR VecSum(vol_vec, &vol);
    CHKERR PetscPrintf(PETSC_COMM_WORLD, "vol  = %3.8e\n", vol);
    if (fabs(vol - 1.0) > 1e-4)
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "Failed to pass test");
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}
