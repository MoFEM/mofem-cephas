/** \file forces_and_sources_testing_users_base.cpp
 * \example forces_and_sources_testing_users_base.cpp
 *
 *  Primarily this is used for testing if the code can handle user base. It is
 *  also, an example of how to build and use user approximation base.
 *  This test calculates CGG base bubble functions and printing them on the
 *  screen. The output of this test is compared with the blessed file.
 *
*/

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#include <MoFEM.hpp>
#include <CGGTonsorialBubbleBase.hpp>

namespace bio = boost::iostreams;
using bio::stream;
using bio::tee_device;

using namespace MoFEM;

static char help[] = "...\n\n";

/**
 * @brief Class used to calculate base functions at integration points
 * 
 */
struct CGGUserPolynomialBase : public BaseFunction {

  CGGUserPolynomialBase() {}
  ~CGGUserPolynomialBase() {}

  /** 
   * @brief Return interface to this class when one ask for for tetrahedron,
   * otherisw return interface class for generic class.
   * 
   * @param uuid Unique id
   * @param iface interface class
   * @return MoFEMErrorCode
   */
   
  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 MoFEM::UnknownInterface **iface) const {
    MoFEMFunctionBegin;
    *iface = NULL;
    if (uuid == IDD_TET_BASE_FUNCTION) {
      *iface = const_cast<CGGUserPolynomialBase *>(this);
      MoFEMFunctionReturnHot(0);
    } else {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "wrong interference");
    }
    CHKERR BaseFunction::query_interface(uuid, iface);
    MoFEMFunctionReturn(0);
  }

  /**
   * @brief Calculate base functions at intergeneration points
   * 
   * @param pts 
   * @param ctx_ptr 
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {
    MoFEMFunctionBeginHot;

    MoFEM::UnknownInterface *iface;
    CHKERR ctx_ptr->query_interface(IDD_TET_BASE_FUNCTION, &iface);
    cTx = reinterpret_cast<EntPolynomialBaseCtx *>(iface);

    int nb_gauss_pts = pts.size2();
    if (!nb_gauss_pts) {
      MoFEMFunctionReturnHot(0);
    }

    if (pts.size1() < 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong dimension of pts, should be at least 3 rows with "
              "coordinates");
    }

    switch (cTx->sPace) {
    case HDIV:
      CHKERR getValueHdivForCGGBubble(pts);
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
    }

    MoFEMFunctionReturnHot(0);
  }

private:
  EntPolynomialBaseCtx *cTx;

  MatrixDouble shapeFun;
  MatrixDouble l2Phi;
  MatrixDouble l2DiffPhi;

  MoFEMErrorCode getValueHdivForCGGBubble(MatrixDouble &pts) {
    MoFEMFunctionBegin;

    const FieldApproximationBase base = cTx->bAse;
    // This should be used only in case USER_BASE is selected
    if(cTx->bAse!=USER_BASE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong base, should be USER_BASE");
    }
    // get access to data structures on element
    DataForcesAndSourcesCore &data = cTx->dAta;
    // Get approximation order on element. Note that bubble functions are only
    // on tetrahedron.
    const int order = data.dataOnEntities[MBTET][0].getDataOrder();
    /// number of integration points
    const int nb_gauss_pts = pts.size2();
    // number of base functions
    const int nb_base_functions = 3 * NBVOLUMETET_CCG_BUBBLE(order);
    // get base functions and set size 
    MatrixDouble &phi = data.dataOnEntities[MBTET][0].getN(USER_BASE);
    phi.resize(nb_gauss_pts, nb_base_functions, false);

    // calculate shape functions, i.e. barycentric coordinates
    shapeFun.resize(nb_gauss_pts, 4, false);
    CHKERR ShapeMBTET(&*shapeFun.data().begin(), &pts(0, 0),
                      &pts(1, 0), &pts(2, 0), nb_gauss_pts);
    // direvatives of shape functions
    double diff_shape_fun[12];
    CHKERR ShapeDiffMBTET(diff_shape_fun);

    // calculate L2 base used to generate bubble HDiv base by CGG recipe
    l2Phi.resize(nb_gauss_pts, NBVOLUMETET_L2(order-2),false);
    l2DiffPhi.resize(nb_gauss_pts, 3 * NBVOLUMETET_L2(order - 2), false);
    CHKERR L2_Ainsworth_ShapeFunctions_MBTET(
        order-2, &shapeFun(0, 0), diff_shape_fun, &l2Phi(0, 0), &l2DiffPhi(0, 0),
        nb_gauss_pts, Legendre_polynomials);

    // finally calculate base functions
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_phi(
        &phi(0, 0), &phi(0, 1), &phi(0, 2));
    CHKERR CGG_BubbleBase_MBTET(order, &shapeFun(0, 0), diff_shape_fun,
                                &l2Phi(0, 0), &l2DiffPhi(0, 0), t_phi,
                                nb_gauss_pts);

    MoFEMFunctionReturn(0);
  }
};

int main(int argc, char *argv[]) {

  // Initialise MoFEM, MPI and petsc 
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    // create moab 
    moab::Core mb_instance;
    // get interface to moab databse
    moab::Interface &moab = mb_instance;

    // get file
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
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "*** ERROR -my_file (MESH FILE NEEDED)");
    }

    // create MoFEM database
    MoFEM::Core core(moab);
    // get interface to moab database
    MoFEM::Interface &m_field = core;

    // load mesh file
    const char *option;
    option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // set bit refinement level
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, BitRefLevel().set(0));

    // Create fields, field "FIELD_CGG" has user base, it means that recipe how
    // to construct approximation is provided by user. Is set that user provided
    // base is in h-div space.
    CHKERR m_field.add_field("FILED_CGG", HDIV, USER_BASE, 1);
    CHKERR m_field.add_field("FILED_RT", HDIV, DEMKOWICZ_JACOBI_BASE, 1);

    // get access to "FIELD_CGG" data structure
    auto field_ptr = m_field.get_field_structure("FILED_CGG");\
    // get table associating number of dofs to entities depending on 
    // approximation order set on those entities.
    auto field_order_table =
        const_cast<Field *>(field_ptr)->getFieldOrderTable();

    // function set zero number of dofs 
    auto get_cgg_bubble_order_zero = [](int p) { return 0; };
    // function set non-zero number of dofs on tetrahedrons
    auto get_cgg_bubble_order_tet = [](int p) {
      return NBVOLUMETET_CCG_BUBBLE(p);
    };
    field_order_table[MBVERTEX] = get_cgg_bubble_order_zero;
    field_order_table[MBEDGE] = get_cgg_bubble_order_zero;
    field_order_table[MBTRI] = get_cgg_bubble_order_zero;
    field_order_table[MBTET] = get_cgg_bubble_order_tet;

    // add finite element
    CHKERR m_field.add_finite_element("FE");

    // define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("FE", "FILED_CGG");
    CHKERR m_field.modify_finite_element_add_field_col("FE", "FILED_CGG");
    CHKERR m_field.modify_finite_element_add_field_data("FE", "FILED_CGG");
    CHKERR m_field.modify_finite_element_add_field_row("FE", "FILED_RT");
    CHKERR m_field.modify_finite_element_add_field_col("FE", "FILED_RT");
    CHKERR m_field.modify_finite_element_add_field_data("FE", "FILED_RT");

    // add problem
    CHKERR m_field.add_problem("PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("PROBLEM", "FE");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("PROBLEM",
                                                    BitRefLevel().set(0));

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FILED_CGG");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FILED_RT");
    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET,
                                                      "FE");

    // set app. order
    int order = 3;
    CHKERR m_field.set_field_order(root_set, MBTRI, "FILED_CGG", order);
    CHKERR m_field.set_field_order(root_set, MBTET, "FILED_CGG", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "FILED_RT", order);
    CHKERR m_field.set_field_order(root_set, MBTET, "FILED_RT", order);

    /****/
    // build database
    // build field
    CHKERR m_field.build_fields();
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(BitRefLevel().set(0));

    // build problem
    CHKERR m_field.getInterface<ProblemsManager>()->buildProblem("PROBLEM",
                                                                 true);
    // dofs partitioning
    CHKERR m_field.getInterface<ProblemsManager>()->partitionSimpleProblem(
        "PROBLEM");
    CHKERR m_field.getInterface<ProblemsManager>()->partitionFiniteElements(
        "PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR m_field.getInterface<ProblemsManager>()->partitionGhostDofs(
        "PROBLEM");

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;

    std::ofstream ofs("forces_and_sources_testing_users_base.txt");
    TeeDevice my_tee(std::cout, ofs);
    TeeStream my_split(my_tee);

    /**
     * Simple user data operator which main purpose is to print values
     * of base functions at intergation points.
     * 
     */
    struct MyOp1 : public VolumeElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &my_split;
      MyOp1(const std::string &row_filed, const std::string &col_field,
            TeeStream &_my_split, char type)
          : VolumeElementForcesAndSourcesCore::UserDataOperator(
                row_filed, col_field, type),
            my_split(_my_split) {
        sYmm = false;
      }

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBeginHot;
        if(data.getIndices().empty()) {
          MoFEMFunctionReturnHot(0);
        }
        my_split << rowFieldName << endl;
        my_split << "side: " << side << " type: " << type << std::endl;
        my_split << data << endl;
        my_split << data.getN() << endl;
        my_split << endl;
        MoFEMFunctionReturnHot(0);
      }

      MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                            EntityType col_type,
                            DataForcesAndSourcesCore::EntData &row_data,
                            DataForcesAndSourcesCore::EntData &col_data) {
        MoFEMFunctionBeginHot;
        if(row_data.getIndices().empty())
          MoFEMFunctionReturnHot(0);
        if(col_data.getIndices().empty())
          MoFEMFunctionReturnHot(0);
        my_split << rowFieldName << " : " << colFieldName << endl;
        my_split << "row side: " << row_side << " row_type: " << row_type
                 << std::endl;
        my_split << "col side: " << col_side << " col_type: " << col_type
                 << std::endl;
        my_split << row_data.getIndices().size() << " : "
                 << col_data.getIndices().size() << endl;
        my_split << endl;
        MoFEMFunctionReturnHot(0);
      }
    };

    // create finite element instance
    VolumeElementForcesAndSourcesCore fe1(m_field);
    // set class needed to cinstruct user approximation base
    fe1.getUserPolynomialBase() =
        boost::shared_ptr<BaseFunction>(new CGGUserPolynomialBase());

    // push user data oprators
    fe1.getOpPtrVector().push_back(
        new MyOp1("FILED_CGG", "FILED_CGG", my_split,
                  ForcesAndSourcesCore::UserDataOperator::OPROW));
    fe1.getOpPtrVector().push_back(
        new MyOp1("FILED_CGG", "FILED_RT", my_split,
                  ForcesAndSourcesCore::UserDataOperator::OPROWCOL));

    // iterate over finite elements, and execute user data operators on each 
    // of them
    CHKERR m_field.loop_finite_elements("PROBLEM", "FE", fe1);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

  return 0;
}