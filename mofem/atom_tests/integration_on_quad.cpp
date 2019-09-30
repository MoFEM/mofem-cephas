/** \file integration_on_quad.cpp
  \example integration_on_quad.cpp
  \brief Testing higher order shape function and 
        their derivatives in integration on quad

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

using namespace MoFEM;
static char help[] = "...\n\n";
static int debug = 1;

using FaceElementUserDataOperator =
    MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator;
using EntData = DataForcesAndSourcesCore::EntData;

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);
  try {
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

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

    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 2, bit_level0);

    // Fields
    CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.add_ents_to_field_by_type(0, MBQUAD, "FIELD1");

    int order = 4;

    CHKERR m_field.set_field_order(0, MBVERTEX, "FIELD1", 1);
    CHKERR m_field.set_field_order(0, MBEDGE, "FIELD1", order);
    CHKERR m_field.set_field_order(0, MBQUAD, "FIELD1", order);

    CHKERR m_field.build_fields();

    // auto set_field = [&](VectorAdaptor &&field_data, double *x, double *y,
    //                      double *z) {
    //   MoFEMFunctionBegin;
    //   // field_data[0] = (-7 * quad(*x) + 7 * cube(*x) + 5 * sqr(*x) - 3 * (*x) - 1) *
    //   //                 (-7 * quad(*y) + 2 * cube(*y) - sqr(*y) + 9 * (*y) - 9);
    //   field_data[0] = sin(M_PI/2.0 * (*x)) * cos(M_PI/2.0 * (*y));
    //   cout << field_data[0];
    //   MoFEMFunctionReturn(0);
    // };

    // FE
    CHKERR m_field.add_finite_element("TEST_FE1");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE1", "FIELD1");

    CHKERR m_field.add_ents_to_finite_element_by_type(0, MBQUAD,
                                                      "TEST_FE1");

    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // //build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM",
                                                     "TEST_FE1");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);

    // build problem
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);
    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    struct CommonData {
     double fieldInteg;
     boost::shared_ptr <VectorDouble> fieldPtr;
     //VectorDouble3 integralFieldGrad;
     CommonData() {
       fieldInteg = 0.0;
       fieldPtr = boost::shared_ptr<VectorDouble>(new VectorDouble());
       //integralFieldGrad.clear();
      }
    };

    struct FieldFunction {
      double operator()(const double x, const double y) const {
        auto sqr = [](double x) { return x * x; };
        auto cube = [](double x) { return x * x * x; };
        auto quad = [](double x) { return x * x * x * x; };
        return (quad(x) * quad(x) + 1) *
               (quad(y) * quad(y) + 1);
        ;
      }
    };

    struct OpCalcIntegral : public FaceElementUserDataOperator {

      typedef boost::function<double(const double, const double)> FieldFunction;

      boost::shared_ptr<CommonData> &commonData;
      FieldFunction fieldFunction;

      OpCalcIntegral(FieldFunction field_function,
                     boost::shared_ptr<CommonData> &common_data)
          : FaceElementUserDataOperator("FIELD1",
                                        FaceElementUserDataOperator::OPROW),
            fieldFunction(field_function), commonData(common_data) {
        doVertices = true;
        doEdges = false;
        doQuads = false;
        doTris = false;
        doTets = false;
        doPrisms = false;
      };

      MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
        MoFEMFunctionBegin;
        if (type != MBVERTEX)
          PetscFunctionReturn(0);

        const int nb_gauss_pts = getGaussPts().size2();
        auto t_coords = getFTensor1CoordsAtGaussPts();

        FTensor::Number<0> NX;
        FTensor::Number<1> NY;

        for (int gg = 0; gg != nb_gauss_pts; gg++) {

          double w = getArea() * getGaussPts()(2, gg);

          commonData->fieldInteg +=
              w * fieldFunction(t_coords(NX), t_coords(NY));

          ++t_coords;
        }

        MoFEMFunctionReturn(0);
      }
    };

    boost::shared_ptr<CommonData> common_data =
        boost::make_shared<CommonData>();

    FaceElementForcesAndSourcesCore fe1(m_field);
    fe1.getOpPtrVector().push_back(
        new OpCalcIntegral(FieldFunction(), common_data));

    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE1", fe1);
    cout.precision(12);
    cout << "Field Integral: " << common_data->fieldInteg << endl;

    //EntityHandle root_set = moab.get_root_set();
    //CHKERR moab.write_file("sphere.vtk", "VTK", "", &root_set, 1);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
